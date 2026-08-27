#!/usr/bin/env bash

set -uo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 RUN_DIRECTORY EXECUTABLE [APPLICATION_ARGUMENT ...]" >&2
  exit 2
fi

script_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
run_directory=$(cd "$1" && pwd)
executable=$2
shift 2
application_arguments=("$@")

ranks=${PRF_BENCH_RANKS:-8}
repetitions=${PRF_BENCH_REPETITIONS:-3}
launcher=${PRF_BENCH_LAUNCHER:-srun}
ompi_info_command=${PRF_BENCH_OMPI_INFO:-ompi_info}
result_directory=${PRF_BENCH_RESULTS:-$run_directory/mpi_prf_benchmark_results}
active_deck=${PRF_BENCH_ACTIVE_DECK:-_data/q2p1_param.dat}
write_deck=${PRF_BENCH_WRITE_DECK:-}
restart_deck=${PRF_BENCH_RESTART_DECK:-}
restart_component=${PRF_BENCH_RESTART_COMPONENT:-}
restart_cb_buffer=${PRF_BENCH_RESTART_CB_BUFFER:-unset}
restart_ds_write=${PRF_BENCH_RESTART_DS_WRITE:-unset}
staged_restart_slot=${PRF_BENCH_STAGED_RESTART_SLOT:-}
success_codes=${PRF_BENCH_SUCCESS_CODES:-0}
print_hints=${PRF_BENCH_PRINT_HINTS:-no}
print_hints_profile=${PRF_BENCH_PRINT_HINTS_PROFILE:-romio_4m_enable}
if [[ $print_hints == 1 ]]; then
  print_hints=yes
fi

read -r -a profiles <<< "${PRF_BENCH_PROFILES:-openmpi_default romio_builtin_defaults romio_4m_enable romio_16m_enable romio_64m_enable romio_16m_disable}"

resolve_run_path() {
  case $1 in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s/%s\n' "$run_directory" "$1" ;;
  esac
}

is_success_code() {
  local actual=$1
  local allowed
  local allowed_codes=${success_codes//,/ }
  for allowed in $allowed_codes; do
    [[ $actual == "$allowed" ]] && return 0
  done
  return 1
}

active_deck_path=$(resolve_run_path "$active_deck")
write_deck_path=
restart_deck_path=
validate_restart=no
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
fi

if [[ -n $staged_restart_slot && $validate_restart != yes ]]; then
  echo "A staged restart slot requires write and restart parameter decks." >&2
  exit 2
fi

if [[ -z ${OMPI_MCA_io:-} ]] && \
    ! command -v "$ompi_info_command" >/dev/null 2>&1; then
  echo "Cannot find ompi_info; set PRF_BENCH_OMPI_INFO or OMPI_MCA_io." >&2
  exit 2
fi

if [[ -n ${OMPI_MCA_io:-} ]]; then
  romio_component=$OMPI_MCA_io
else
  romio_component=$("$ompi_info_command" --param io all 2>/dev/null |
    sed -n 's/.*MCA io: \(romio[[:alnum:]_]*\).*/\1/p' | head -1)
fi
if [[ -z $romio_component ]]; then
  echo "No OpenMPI ROMIO component detected; set OMPI_MCA_io explicitly." >&2
  exit 2
fi

mkdir -p "$result_directory"
result_directory=$(cd "$result_directory" && pwd)
hints_directory="$result_directory/romio_hints"
mkdir -p "$hints_directory"
summary="$result_directory/summary.tsv"
printf 'ranks\tprofile\trepetition\tcomponent\tromio_hints\thints_source\tcb_buffer_size\tromio_ds_write\tcheckpoint_s\tpayload_bytes\tthroughput_mib_s\trun_elapsed_s\twrite_exit_code\twrite_status\treload_elapsed_s\treload_exit_code\troundtrip_status\n' > "$summary"
profile_settings="$result_directory/profile_settings.tsv"
printf 'profile\tcomponent\tromio_hints\thints_source\tcb_buffer_size\tromio_ds_write\n' > "$profile_settings"

cleanup() {
  local exit_code=$?
  if [[ $validate_restart == yes && -f $write_deck_path ]]; then
    cp "$write_deck_path" "$active_deck_path"
  fi
  if [[ -s $summary ]]; then
    python3 "$script_directory/summarize_romio_matrix.py" "$summary" || true
  fi
  exit "$exit_code"
}
trap cleanup EXIT

{
  printf 'timestamp_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'host=%s\n' "$(hostname)"
  printf 'ranks=%s\n' "$ranks"
  printf 'repetitions=%s\n' "$repetitions"
  printf 'launcher=%s\n' "$launcher"
  printf 'executable=%s\n' "$executable"
  printf 'romio_component=%s\n' "$romio_component"
  printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-unset}"
  if command -v "$ompi_info_command" >/dev/null 2>&1; then
    "$ompi_info_command" --version 2>/dev/null | head -1 || true
    "$ompi_info_command" --param io all 2>/dev/null || true
  fi
} > "$result_directory/environment.txt"

run_application() {
  local selected_component=$1
  local selected_hints=$2
  local selected_print_hints=${3:-no}
  local print_hints_env=(-u ROMIO_PRINT_HINTS)

  if [[ $selected_print_hints == yes ]]; then
    print_hints_env=(ROMIO_PRINT_HINTS=1)
  fi

  if [[ $selected_component == unset ]]; then
    env -u OMPI_MCA_io -u ROMIO_HINTS -u ROMIO_CB_BUFFER_SIZE \
      -u ROMIO_DS_WRITE "${print_hints_env[@]}" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  elif [[ $selected_hints == unset ]]; then
    env -u ROMIO_HINTS -u ROMIO_CB_BUFFER_SIZE -u ROMIO_DS_WRITE \
      "${print_hints_env[@]}" \
      OMPI_MCA_io="$selected_component" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  else
    env -u ROMIO_CB_BUFFER_SIZE -u ROMIO_DS_WRITE \
      "${print_hints_env[@]}" OMPI_MCA_io="$selected_component" \
      ROMIO_HINTS="$selected_hints" \
      "$launcher" --ntasks="$ranks" "$executable" \
      "${application_arguments[@]}"
  fi
}

write_romio_hints() {
  local hints_file=$1
  local selected_cb_buffer=$2
  local selected_ds_write=$3

  : > "$hints_file"
  if [[ $selected_cb_buffer != unset ]]; then
    printf 'cb_buffer_size %s\n' "$selected_cb_buffer" >> "$hints_file"
  fi
  if [[ $selected_ds_write != unset ]]; then
    printf 'romio_ds_write %s\n' "$selected_ds_write" >> "$hints_file"
  fi
}

verify_printed_hint() {
  local log_file=$1
  local key=$2
  local expected_value=$3

  awk -v key="$key" -v expected="$expected_value" '
    index($0, key) && (expected == "unset" || index($0, expected)) { found=1 }
    END { exit !found }
  ' "$log_file"
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

hint_probe_done=no
for profile in "${profiles[@]}"; do
  component=$romio_component
  profile_hints=unset
  hints_source=none
  cb_buffer_size=unset
  romio_ds_write=unset
  case "$profile" in
    openmpi_default)
      component=unset
      hints_source=openmpi_default
      ;;
    romio_builtin_defaults)
      hints_source=romio_builtin_defaults
      ;;
    romio_4m_enable)
      cb_buffer_size=4194304
      romio_ds_write=enable
      ;;
    romio_16m_enable)
      cb_buffer_size=16777216
      romio_ds_write=enable
      ;;
    romio_64m_enable)
      cb_buffer_size=67108864
      romio_ds_write=enable
      ;;
    romio_16m_disable)
      cb_buffer_size=16777216
      romio_ds_write=disable
      ;;
    *)
      echo "Unknown profile: $profile" >&2
      exit 2
      ;;
  esac
  if [[ $component != unset && $cb_buffer_size != unset ]] || \
      [[ $component != unset && $romio_ds_write != unset ]]; then
    profile_hints="$hints_directory/${profile}.conf"
    write_romio_hints "$profile_hints" "$cb_buffer_size" "$romio_ds_write"
    hints_source=ROMIO_HINTS
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$profile" "$component" \
    "$profile_hints" "$hints_source" "$cb_buffer_size" "$romio_ds_write" \
    >> "$profile_settings"

  if [[ $print_hints == yes && $profile == "$print_hints_profile" ]]; then
    if [[ $component == unset ]]; then
      echo "ROMIO hint verification requires a ROMIO profile: $profile" >&2
      exit 2
    fi
    if [[ $validate_restart == yes ]]; then
      cp "$write_deck_path" "$active_deck_path"
    fi
    verify_log="$result_directory/${profile}.romio_print_hints.log"
    echo "Running untimed ROMIO_PRINT_HINTS verification for $profile." >&2
    (
      cd "$run_directory" || exit 2
      run_application "$component" "$profile_hints" yes
    ) > "$verify_log" 2>&1
    verify_exit_code=$?
    if ! is_success_code "$verify_exit_code"; then
      echo "ROMIO hint verification failed; see $verify_log" >&2
      exit "$verify_exit_code"
    fi
    if ! verify_printed_hint "$verify_log" cb_buffer_size \
        "$cb_buffer_size" || \
        ! verify_printed_hint "$verify_log" romio_ds_write \
        "$romio_ds_write"; then
      echo "ROMIO did not print the requested hints; see $verify_log" >&2
      exit 1
    fi
    hint_probe_done=yes
  fi

  for ((repetition=1; repetition<=repetitions; repetition++)); do
    log_file="$result_directory/${profile}_np${ranks}_r${repetition}.log"
    restart_log_file="$result_directory/${profile}_np${ranks}_r${repetition}.restart.log"
    checkpoint_archive="$result_directory/checkpoints/${profile}_np${ranks}_r${repetition}"
    if [[ $validate_restart == yes ]]; then
      cp "$write_deck_path" "$active_deck_path"
    fi

    start_ns=$(date +%s%N)
    (
      cd "$run_directory" || exit 2
      run_application "$component" "$profile_hints"
    ) > "$log_file" 2>&1
    write_exit_code=$?
    end_ns=$(date +%s%N)
    run_elapsed=$(awk -v start="$start_ns" -v end="$end_ns" \
      'BEGIN { printf "%.6f", (end-start)/1000000000 }')
    if is_success_code "$write_exit_code"; then
      write_status=ok
    else
      write_status=failed
    fi

    checkpoint_slot=
    checkpoint_elapsed=not_available
    checkpoint_directory=
    payload_bytes=0
    throughput=not_available
    reload_elapsed=not_run
    reload_exit_code=not_run
    roundtrip_status=not_run
    failure_message=
    failure_exit_code=0

    if [[ $write_status != ok ]]; then
      failure_message="$profile repetition $repetition solver write failed"
      failure_message+=" with exit code $write_exit_code; see $log_file"
      failure_exit_code=$write_exit_code
    else
      parsed_checkpoint_slot=$(awk '/MPI-PRF checkpoint slot/ {
        for (i=1; i<=NF; i++) if ($i == "slot") slot=$(i+1)
      } END { if (slot != "") print slot }' "$log_file")
      parsed_checkpoint_elapsed=$(awk '/MPI-PRF checkpoint slot/ {
        for (i=1; i<=NF; i++) if ($i == "in") elapsed=$(i+1)
      } END { if (elapsed != "") print elapsed }' "$log_file")
      if [[ -z $parsed_checkpoint_slot || -z $parsed_checkpoint_elapsed ]]; then
        failure_message="$profile repetition $repetition produced no timed checkpoint."
        failure_exit_code=1
      else
        checkpoint_slot=$parsed_checkpoint_slot
        checkpoint_elapsed=$parsed_checkpoint_elapsed
        checkpoint_directory="$run_directory/_dump/$checkpoint_slot"
        if [[ ! -d $checkpoint_directory ]]; then
          failure_message="Timed checkpoint directory does not exist: $run_directory/_dump/$checkpoint_slot"
          failure_exit_code=1
        fi
      fi
      if [[ -z $failure_message ]]; then
        payload_bytes=$(find "$checkpoint_directory" -maxdepth 1 -type f \
          \( -name '*_key*.prf' -o -name '*_comp*.prf' \) -printf '%s\n' |
          awk '{total += $1} END {print total+0}')
        if [[ $payload_bytes -eq 0 ]]; then
          failure_message="Timed checkpoint contains no PRF payload files."
          failure_exit_code=1
        else
          throughput=$(awk -v bytes="$payload_bytes" -v seconds="$checkpoint_elapsed" \
            'BEGIN { printf "%.6f", bytes/1048576/seconds }')
        fi
      fi
    fi

    if [[ -n $failure_message && $write_status == ok ]]; then
      write_status=checkpoint_failed
    fi

    if [[ -z $failure_message && $validate_restart == yes ]]; then
      mkdir -p "$(dirname "$checkpoint_archive")"
      cp -a "$checkpoint_directory" "$checkpoint_archive"
      if [[ -n $staged_restart_slot ]]; then
        staged_directory="$run_directory/_dump/$staged_restart_slot"
        if [[ -e $staged_directory ]]; then
          echo "Staged restart slot already exists: $staged_directory" >&2
          exit 1
        fi
        mv "$checkpoint_directory" "$staged_directory"
      fi

      cp "$restart_deck_path" "$active_deck_path"
      reload_start_ns=$(date +%s%N)
      (
        cd "$run_directory" || exit 2
        if [[ -z $restart_component ]]; then
          run_application "$component" "$profile_hints"
        else
          restart_hints=unset
          if [[ $restart_cb_buffer != unset || $restart_ds_write != unset ]]; then
            restart_hints="$hints_directory/${profile}.restart.conf"
            write_romio_hints "$restart_hints" "$restart_cb_buffer" \
              "$restart_ds_write"
          fi
          run_application "$restart_component" "$restart_hints"
        fi
      ) > "$restart_log_file" 2>&1
      reload_exit_code=$?
      reload_end_ns=$(date +%s%N)
      reload_elapsed=$(awk -v start="$reload_start_ns" -v end="$reload_end_ns" \
        'BEGIN { printf "%.6f", (end-start)/1000000000 }')

      if is_success_code "$reload_exit_code"; then
        if [[ -n $staged_restart_slot ]]; then
          if compare_checkpoint_payloads "$checkpoint_archive" \
              "$checkpoint_directory"; then
            roundtrip_status=bitwise_equal
            rm -r "$staged_directory"
          else
            roundtrip_status=mismatch
          fi
        else
          roundtrip_status=reload_ok
        fi
      else
        roundtrip_status=reload_failed
      fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ranks" "$profile" "$repetition" "$component" "$profile_hints" \
      "$hints_source" "$cb_buffer_size" "$romio_ds_write" \
      "$checkpoint_elapsed" "$payload_bytes" "$throughput" "$run_elapsed" \
      "$write_exit_code" "$write_status" "$reload_elapsed" \
      "$reload_exit_code" "$roundtrip_status" >> "$summary"

    grep 'MPI-PRF checkpoint slot' "$log_file" | tail -1 || true
    if [[ -n $failure_message ]]; then
      echo "$failure_message" >&2
      exit "$failure_exit_code"
    fi
    if [[ $validate_restart == yes ]] && \
        ! is_success_code "$reload_exit_code"; then
      echo "$profile repetition $repetition failed restart validation;" \
        "see $restart_log_file" >&2
      exit "$reload_exit_code"
    fi
    if [[ -n $staged_restart_slot && \
        $roundtrip_status != bitwise_equal ]]; then
      echo "$profile repetition $repetition failed bitwise round trip." >&2
      exit 1
    fi
  done
done

if [[ $print_hints == yes && $hint_probe_done != yes ]]; then
  echo "ROMIO hint verification profile is not in PRF_BENCH_PROFILES:" \
    "$print_hints_profile" >&2
  exit 2
fi

echo "Wrote $summary"
