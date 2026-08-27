#!/usr/bin/env bash

set -u

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 RUN_DIRECTORY EXECUTABLE [APPLICATION_ARGUMENT ...]" >&2
  exit 2
fi

run_directory=$1
executable=$2
shift 2
application_arguments=("$@")

ranks=${PRF_BENCH_RANKS:-8}
repetitions=${PRF_BENCH_REPETITIONS:-3}
launcher=${PRF_BENCH_LAUNCHER:-srun}
result_directory=${PRF_BENCH_RESULTS:-mpi_prf_benchmark_results}
active_deck=${PRF_BENCH_ACTIVE_DECK:-_data/q2p1_param.dat}
write_deck=${PRF_BENCH_WRITE_DECK:-}
restart_deck=${PRF_BENCH_RESTART_DECK:-}
restart_component=${PRF_BENCH_RESTART_COMPONENT:-}
restart_cb_buffer=${PRF_BENCH_RESTART_CB_BUFFER:-unset}
restart_ds_write=${PRF_BENCH_RESTART_DS_WRITE:-unset}
write_slot=${PRF_BENCH_WRITE_SLOT:-1}
staged_restart_slot=${PRF_BENCH_STAGED_RESTART_SLOT:-}

if [[ ! -d $run_directory ]]; then
  echo "Run directory does not exist: $run_directory" >&2
  exit 2
fi

resolve_run_path() {
  case $1 in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$run_directory" "$1" ;;
  esac
}

active_deck_path=$(resolve_run_path "$active_deck")
if [[ -n $write_deck || -n $restart_deck ]]; then
  if [[ -z $write_deck || -z $restart_deck ]]; then
    echo "Set both PRF_BENCH_WRITE_DECK and PRF_BENCH_RESTART_DECK." >&2
    exit 2
  fi
  write_deck_path=$(resolve_run_path "$write_deck")
  restart_deck_path=$(resolve_run_path "$restart_deck")
  if [[ ! -f $write_deck_path || ! -f $restart_deck_path ]]; then
    echo "Write or restart parameter template does not exist." >&2
    exit 2
  fi
  validate_restart=yes
  trap 'cp "$write_deck_path" "$active_deck_path"' EXIT
else
  validate_restart=no
fi

if [[ -n $staged_restart_slot && $validate_restart != yes ]]; then
  echo "A staged restart slot requires write and restart parameter decks." >&2
  exit 2
fi

if [[ -n ${OMPI_MCA_io:-} ]]; then
  romio_component=$OMPI_MCA_io
else
  romio_component=$(ompi_info --param io all 2>/dev/null |
    sed -n 's/.*MCA io: \(romio[[:alnum:]_]*\).*/\1/p' | head -1)
fi
if [[ -z $romio_component ]]; then
  echo "No OpenMPI ROMIO component detected; set OMPI_MCA_io explicitly." >&2
  exit 2
fi

mkdir -p "$result_directory"
summary="$result_directory/summary.tsv"
printf 'ranks\tprofile\trepetition\tcomponent\tcb_buffer\tds_write\telapsed_s\texit_code\tcheckpoint_s\tpayload_bytes\treload_elapsed_s\treload_exit_code\troundtrip_status\n' > "$summary"

profiles=(
  openmpi_default
  romio_builtin_defaults
  romio_4m_enable
  romio_16m_enable
  romio_64m_enable
  romio_16m_disable
)

run_profile() {
  if [[ $component == unset ]]; then
    env -u OMPI_MCA_io -u ROMIO_CB_BUFFER_SIZE -u ROMIO_DS_WRITE \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  elif [[ $cb_buffer == unset ]]; then
    env -u ROMIO_CB_BUFFER_SIZE -u ROMIO_DS_WRITE \
      OMPI_MCA_io="$component" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  else
    env OMPI_MCA_io="$component" \
      ROMIO_CB_BUFFER_SIZE="$cb_buffer" ROMIO_DS_WRITE="$ds_write" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  fi
}

run_restart() {
  if [[ -z $restart_component ]]; then
    run_profile
  elif [[ $restart_cb_buffer == unset ]]; then
    env -u ROMIO_CB_BUFFER_SIZE -u ROMIO_DS_WRITE \
      OMPI_MCA_io="$restart_component" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  else
    env OMPI_MCA_io="$restart_component" \
      ROMIO_CB_BUFFER_SIZE="$restart_cb_buffer" \
      ROMIO_DS_WRITE="$restart_ds_write" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  fi
}

compare_checkpoint_payloads() {
  local original=$1
  local regenerated=$2
  local original_files=()
  local regenerated_files=()
  local index name

  shopt -s nullglob
  original_files=("$original"/*_key*.prf "$original"/*_comp*.prf)
  regenerated_files=("$regenerated"/*_key*.prf \
    "$regenerated"/*_comp*.prf)
  shopt -u nullglob
  if [[ ${#original_files[@]} -eq 0 || \
      ${#original_files[@]} -ne ${#regenerated_files[@]} ]]; then
    return 1
  fi
  for ((index=0; index<${#original_files[@]}; index++)); do
    name=$(basename "${original_files[index]}")
    if [[ $name != $(basename "${regenerated_files[index]}") ]] || \
        ! cmp -s "${original_files[index]}" \
          "${regenerated_files[index]}"; then
      return 1
    fi
  done
}

for profile in "${profiles[@]}"; do
  component=$romio_component
  cb_buffer=unset
  ds_write=unset
  case "$profile" in
    openmpi_default)
      component=unset
      ;;
    romio_builtin_defaults)
      ;;
    romio_4m_enable)
      cb_buffer=4194304
      ds_write=enable
      ;;
    romio_16m_enable)
      cb_buffer=16777216
      ds_write=enable
      ;;
    romio_64m_enable)
      cb_buffer=67108864
      ds_write=enable
      ;;
    romio_16m_disable)
      cb_buffer=16777216
      ds_write=disable
      ;;
  esac

  for ((repetition=1; repetition<=repetitions; repetition++)); do
    log_file="$result_directory/${profile}_np${ranks}_r${repetition}.log"
    restart_log_file="$result_directory/${profile}_np${ranks}_r${repetition}.restart.log"
    if [[ $validate_restart == yes ]]; then
      cp "$write_deck_path" "$active_deck_path"
    fi
    start_ns=$(date +%s%N)
    (
      cd "$run_directory" || exit 2
      run_profile
    ) > "$log_file" 2>&1
    exit_code=$?
    end_ns=$(date +%s%N)
    elapsed=$(awk -v start="$start_ns" -v end="$end_ns" \
      'BEGIN { printf "%.6f", (end-start)/1000000000 }')
    checkpoint_elapsed=$(awk '/MPI-PRF checkpoint slot/ {
      for (i=1; i<=NF; i++) if ($i == "in") elapsed=$(i+1)
    } END { if (elapsed != "") print elapsed }' "$log_file")
    payload_bytes=$(awk '/MPI-PRF checkpoint slot/ {
      for (i=1; i<=NF; i++) if ($i == "payload") bytes=$(i-1)
    } END { gsub(/[()]/,"",bytes); if (bytes != "") print bytes }' \
      "$log_file")
    checkpoint_elapsed=${checkpoint_elapsed:-missing}
    payload_bytes=${payload_bytes:-missing}
    reload_elapsed=not_run
    reload_exit_code=not_run
    roundtrip_status=not_run

    if [[ $exit_code -eq 0 && $checkpoint_elapsed != missing && \
        $validate_restart == yes ]]; then
      if [[ -n $staged_restart_slot ]]; then
        original_slot="$run_directory/_dump/$write_slot"
        staged_slot="$run_directory/_dump/$staged_restart_slot"
        if [[ ! -d $original_slot || -e $staged_slot ]]; then
          echo "Cannot stage checkpoint $original_slot as $staged_slot." >&2
          exit 1
        fi
        mv "$original_slot" "$staged_slot"
      fi
      cp "$restart_deck_path" "$active_deck_path"
      reload_start_ns=$(date +%s%N)
      (
        cd "$run_directory" || exit 2
        run_restart
      ) > "$restart_log_file" 2>&1
      reload_exit_code=$?
      reload_end_ns=$(date +%s%N)
      reload_elapsed=$(awk -v start="$reload_start_ns" -v end="$reload_end_ns" \
        'BEGIN { printf "%.6f", (end-start)/1000000000 }')
      if [[ $reload_exit_code -eq 0 && -n $staged_restart_slot ]]; then
        regenerated_slot="$run_directory/_dump/$write_slot"
        if compare_checkpoint_payloads "$staged_slot" "$regenerated_slot"; then
          roundtrip_status=bitwise_equal
          checkpoint_archive="$run_directory/$result_directory/checkpoints"
          mkdir -p "$checkpoint_archive"
          mv "$staged_slot" \
            "$checkpoint_archive/${profile}_np${ranks}_r${repetition}"
        else
          roundtrip_status=mismatch
        fi
      elif [[ $reload_exit_code -eq 0 ]]; then
        roundtrip_status=reload_ok
      else
        roundtrip_status=reload_failed
      fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ranks" "$profile" "$repetition" "$component" "$cb_buffer" \
      "$ds_write" "$elapsed" "$exit_code" "$checkpoint_elapsed" \
      "$payload_bytes" "$reload_elapsed" "$reload_exit_code" \
      "$roundtrip_status" >> "$summary"
    grep 'MPI-PRF checkpoint slot' "$log_file" || true
    if [[ $exit_code -ne 0 ]]; then
      echo "$profile repetition $repetition failed; see $log_file" >&2
      exit "$exit_code"
    fi
    if [[ $checkpoint_elapsed == missing ]]; then
      echo "$profile repetition $repetition produced no timed checkpoint." >&2
      exit 1
    fi
    if [[ $validate_restart == yes && $reload_exit_code -ne 0 ]]; then
      echo "$profile repetition $repetition failed restart validation;" \
        "see $restart_log_file" >&2
      exit "$reload_exit_code"
    fi
    if [[ -n $staged_restart_slot && $roundtrip_status != bitwise_equal ]]; then
      echo "$profile repetition $repetition failed bitwise round-trip;" \
        "preserved slots $staged_restart_slot and $write_slot" >&2
      exit 1
    fi
  done
done

echo "Wrote $summary"
