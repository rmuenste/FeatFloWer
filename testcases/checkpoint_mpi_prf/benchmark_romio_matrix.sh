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

if [[ ! -d $run_directory ]]; then
  echo "Run directory does not exist: $run_directory" >&2
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
printf 'ranks\tprofile\trepetition\tcomponent\tcb_buffer\tds_write\telapsed_s\texit_code\n' > "$summary"

profiles=(
  openmpi_default
  romio_builtin_defaults
  romio_4m_enable
  romio_16m_enable
  romio_64m_enable
  romio_16m_disable
)

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
    start_ns=$(date +%s%N)
    (
      cd "$run_directory" || exit 2
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
    ) > "$log_file" 2>&1
    exit_code=$?
    end_ns=$(date +%s%N)
    elapsed=$(awk -v start="$start_ns" -v end="$end_ns" \
      'BEGIN { printf "%.6f", (end-start)/1000000000 }')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ranks" "$profile" "$repetition" "$component" "$cb_buffer" \
      "$ds_write" "$elapsed" "$exit_code" >> "$summary"
    grep 'MPI-PRF checkpoint slot' "$log_file" || true
    if [[ $exit_code -ne 0 ]]; then
      echo "$profile repetition $repetition failed; see $log_file" >&2
    fi
  done
done

echo "Wrote $summary"
