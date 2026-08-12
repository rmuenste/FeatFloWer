#!/bin/bash
#SBATCH --job-name=ff-hypre-sse
#SBATCH --partition=short
#SBATCH --nodelist=tardis
#SBATCH --constraint=ampere
#SBATCH --gres=gpu:ampere:1
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=384G
#SBATCH --time=02:00:00
#SBATCH --output=stage1-work/slurm-%x-%j.out

set -euo pipefail

module purge
module load gcc/13.2.0
module load openmpi/options/interface/ethernet openmpi/options/cuda/no openmpi/4.1.6
module load cuda/12.4

repo_root="${FF_SSE_ROOT:-${SLURM_SUBMIT_DIR}}"
cpu_app="${FF_SSE_CPU_APP:-${repo_root}/stage1-work/build-matched-cpu/applications/q2p1_sse}"
gpu_app="${FF_SSE_GPU_APP:-${repo_root}/stage1-work/build-gpu/applications/q2p1_sse}"
steps="${FF_SSE_STEPS:-12}"
repeats="${FF_SSE_REPEATS:-1}"
warmup_calls="${FF_SSE_WARMUP_CALLS:-0}"
backends="${FF_SSE_BACKENDS:-cpu gpu}"
ranks_list="${FF_SSE_RANKS:-32}"
solver_types="${FF_SSE_SOLVER_TYPES:-7 8}"
mesh_level="${FF_SSE_MESH_LEVEL:-2}"
pressure_min_level="${FF_SSE_PRESSURE_MIN_LEVEL:-1}"
mesh_quality="${FF_SSE_MESH_QUALITY:-coarse}"
use_mps="${FF_SSE_USE_MPS:-1}"
run_id="${FF_SSE_RUN_ID:-sse-${SLURM_JOB_ID}}"
result_root="${repo_root}/stage1-work/results/${run_id}"
project_folder="_ianus/TSE/Conv"

for app_dir in "${cpu_app}" "${gpu_app}"; do
  for required in q2p1_sse s3d_mesher libmetis.so _data/MG.dat \
    _data_BU/q2p1_paramV_BU_test.dat "${project_folder}/setup.e3d"; do
    if [[ ! -e "${app_dir}/${required}" ]]; then
      echo "Missing SSE runtime asset: ${app_dir}/${required}" >&2
      exit 2
    fi
  done
done

mkdir -p "${result_root}/runs" "${result_root}/metadata"

{
  echo "run_id=${run_id}"
  echo "host=$(hostname)"
  echo "job_id=${SLURM_JOB_ID}"
  echo "steps=${steps}"
  echo "repeats=${repeats}"
  echo "warmup_calls=${warmup_calls}"
  echo "backends=${backends}"
  echo "ranks=${ranks_list}"
  echo "solver_types=${solver_types}"
  echo "mesh_level=${mesh_level}"
  echo "pressure_min_level=${pressure_min_level}"
  echo "mesh_quality=${mesh_quality}"
  echo "use_mps=${use_mps}"
  echo "project_folder=${project_folder}"
  echo "cpu_app=${cpu_app}"
  echo "gpu_app=${gpu_app}"
  echo "git_commit=$(git -C "${repo_root}" rev-parse HEAD)"
  echo "cuda=$(/sfw/cuda/12.4/bin/nvcc --version | tail -1)"
  echo "mpi=$(mpirun --version | head -1)"
} > "${result_root}/metadata/run.env"
nvidia-smi -q > "${result_root}/metadata/nvidia-smi.txt"

export CUDA_VISIBLE_DEVICES=0
export PYTHONUNBUFFERED=1
export OMPI_MCA_io=romio321
export ROMIO_CB_BUFFER_SIZE=16777216
export ROMIO_DS_WRITE=enable

if [[ "${use_mps}" == 1 ]]; then
  mps_pipe="${TMPDIR:-/tmp}/ff-sse-mps-${SLURM_JOB_ID}"
  mps_log="${TMPDIR:-/tmp}/ff-sse-mps-log-${SLURM_JOB_ID}"
  mkdir -p "${mps_pipe}" "${mps_log}"
  export CUDA_MPS_PIPE_DIRECTORY="${mps_pipe}"
  export CUDA_MPS_LOG_DIRECTORY="${mps_log}"
  nvidia-cuda-mps-control -d
fi
nvidia-smi dmon -s pucvmet -d 1 \
  > "${result_root}/metadata/nvidia-smi-dmon.log" 2>&1 &
dmon_pid=$!

cleanup() {
  kill "${dmon_pid}" 2>/dev/null || true
  if [[ "${use_mps}" == 1 ]]; then
    echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

stage_case() {
  local backend=$1
  local app_dir=$2
  local ranks=$3
  local solver_type=$4
  local repetition=$5
  local case_dir=$6
  local parameter_file

  mkdir -p "${case_dir}"
  cp -a "${app_dir}/." "${case_dir}/"
  sed -i "s|^MeshQuality=.*|MeshQuality=${mesh_quality}|" \
    "${case_dir}/${project_folder}/setup.e3d"
  parameter_file="${case_dir}/_data_BU/q2p1_paramV_BU_test.dat"
  sed -i \
    -e "s|^SimPar@MaxMeshLevel.*|SimPar@MaxMeshLevel = ${mesh_level}|" \
    -e "s|^SimPar@MaxNumStep.*|SimPar@MaxNumStep = ${steps}|" \
    -e 's|^SimPar@BackUpFreq.*|SimPar@BackUpFreq = 1000000|' \
    -e 's|^SimPar@OutputFreq.*|SimPar@OutputFreq = 1000000d0|' \
    -e "s|^Pres@MGMinLev.*|Pres@MGMinLev = ${pressure_min_level}|" \
    -e "s|^Pres@MGMedLev.*|Pres@MGMedLev = ${pressure_min_level}|" \
    -e "s|^Pres@MGCrsSolverType.*|Pres@MGCrsSolverType = ${solver_type}|" \
    "${parameter_file}"

  {
    echo "backend=${backend}"
    echo "ranks=${ranks}"
    echo "workers=$((ranks - 1))"
    echo "solver_type=${solver_type}"
    echo "repetition=${repetition}"
    echo "mesh_level=${mesh_level}"
    echo "pressure_min_level=${pressure_min_level}"
    echo "mesh_quality=${mesh_quality}"
    echo "steps=${steps}"
  } > "${case_dir}/case.env"
}

for backend in ${backends}; do
  if [[ "${backend}" == cpu ]]; then
    app_dir="${cpu_app}"
  else
    app_dir="${gpu_app}"
  fi

  for ranks in ${ranks_list}; do
    for solver_type in ${solver_types}; do
      for repetition in $(seq 1 "${repeats}"); do
        name="${backend}-np${ranks}-type${solver_type}-rep${repetition}"
        case_dir="${result_root}/runs/${name}"
        stage_case "${backend}" "${app_dir}" "${ranks}" \
          "${solver_type}" "${repetition}" "${case_dir}"
        set +e
        (
          cd "${case_dir}"
          export LD_LIBRARY_PATH="${case_dir}:${LD_LIBRARY_PATH:-}"
          /usr/bin/time -p python3 ./e3d_start.py -n "${ranks}" \
            -f "${project_folder}" -a 0 --short-test \
            --max-meshlevel "${mesh_level}"
        ) > "${case_dir}/run.log" 2>&1
        run_exit_code=$?
        set -e
        echo "run_exit_code=${run_exit_code}" >> "${case_dir}/case.env"
        grep -E 'HYPRE_TIMING|stopped|error|failed|abnormal|CurrentStatus' \
          "${case_dir}/run.log" > "${case_dir}/summary.log" || true
      done
    done
  done
done

python3 "${repo_root}/tools/cluster_scripts/summarize_hypre_gpu_stage1.py" \
  "${result_root}" --warmup-calls "${warmup_calls}" \
  --output "${result_root}/summary.tsv"

echo "Stage 1 SSE results: ${result_root}"
