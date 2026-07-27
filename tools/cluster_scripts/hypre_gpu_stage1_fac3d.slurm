#!/bin/bash
#SBATCH --job-name=ff-hypre-fac3d
#SBATCH --partition=short
#SBATCH --nodelist=tardis
#SBATCH --constraint=ampere
#SBATCH --gres=gpu:ampere:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
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

repo_root="${FF_FAC3D_ROOT:-${SLURM_SUBMIT_DIR}}"
cpu_app="${FF_FAC3D_CPU_APP:-}"
gpu_app="${FF_FAC3D_GPU_APP:-}"
if [[ -z "${cpu_app}" ]]; then
  cpu_app="${repo_root}/stage1-work/build-matched-cpu/applications/q2p1_fac3d"
fi
if [[ -z "${gpu_app}" ]]; then
  gpu_app="${repo_root}/stage1-work/build-gpu/applications/q2p1_fac3d"
fi
steps="${FF_FAC3D_STEPS:-1}"
repeats="${FF_FAC3D_REPEATS:-1}"
warmup_calls="${FF_FAC3D_WARMUP_CALLS:-0}"
mesh_level="${FF_FAC3D_MESH_LEVEL:-5}"
backends="${FF_FAC3D_BACKENDS:-cpu gpu}"
ranks_list="${FF_FAC3D_RANKS:-4}"
solver_types="${FF_FAC3D_SOLVER_TYPES:-7 8}"
use_mps="${FF_FAC3D_USE_MPS:-1}"
keep_output="${FF_FAC3D_KEEP_OUTPUT:-1}"
run_id="${FF_FAC3D_RUN_ID:-fac3d-${SLURM_JOB_ID}}"
result_root="${repo_root}/stage1-work/results/${run_id}"
partition_root="${cpu_app}/_mesh"
project_file="_adc/3D_FAC_FBM/file.prj"

if [[ ! -x "${cpu_app}/q2p1_fac3d" ]]; then
  echo "Missing CPU executable: ${cpu_app}/q2p1_fac3d" >&2
  exit 2
fi
if [[ ! -x "${gpu_app}/q2p1_fac3d" ]]; then
  echo "Missing GPU executable: ${gpu_app}/q2p1_fac3d" >&2
  exit 2
fi
if [[ ! -f "${cpu_app}/${project_file}" ]]; then
  echo "Missing 3D cylinder project: ${cpu_app}/${project_file}" >&2
  exit 2
fi

mkdir -p "${result_root}/runs" "${result_root}/metadata"

{
  echo "run_id=${run_id}"
  echo "host=$(hostname)"
  echo "job_id=${SLURM_JOB_ID}"
  echo "mesh_level=${mesh_level}"
  echo "steps=${steps}"
  echo "repeats=${repeats}"
  echo "warmup_calls=${warmup_calls}"
  echo "backends=${backends}"
  echo "ranks=${ranks_list}"
  echo "solver_types=${solver_types}"
  echo "use_mps=${use_mps}"
  echo "keep_output=${keep_output}"
  echo "project_file=${project_file}"
  echo "git_commit=$(git -C "${repo_root}" rev-parse HEAD)"
  echo "cuda=$(/sfw/cuda/12.4/bin/nvcc --version | tail -1)"
  echo "mpi=$(mpirun --version | head -1)"
} > "${result_root}/metadata/run.env"
nvidia-smi -q > "${result_root}/metadata/nvidia-smi.txt"

for ranks in ${ranks_list}; do
  workers=$((ranks - 1))
  mesh_folder="FAC3D_np${ranks}"
  mesh_file="${partition_root}/${mesh_folder}/sub0001/GRID0001.tri"
  if [[ ! -f "${mesh_file}" ]]; then
    (
      cd "${cpu_app}"
      LD_LIBRARY_PATH="${cpu_app}:${LD_LIBRARY_PATH:-}" \
        python3 "${repo_root}/tools/PyPartitioner.py" \
        "${workers}" 1 1 "${mesh_folder}" "${project_file}"
    ) > "${result_root}/metadata/partition-np${ranks}.log" 2>&1
  fi
  if [[ ! -f "${mesh_file}" ]]; then
    echo "Partitioning did not create ${mesh_file}" >&2
    exit 3
  fi
done

export CUDA_VISIBLE_DEVICES=0

if [[ "${use_mps}" == 1 ]]; then
  mps_pipe="${TMPDIR:-/tmp}/ff-fac3d-mps-${SLURM_JOB_ID}"
  mps_log="${TMPDIR:-/tmp}/ff-fac3d-mps-log-${SLURM_JOB_ID}"
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
  local mesh_folder="FAC3D_np${ranks}"

  mkdir -p "${case_dir}/_data" "${case_dir}/_dump" "${case_dir}/_gmv" \
    "${case_dir}/_mesh" "${case_dir}/_ns" "${case_dir}/_vtk" \
    "${case_dir}/solution" "${case_dir}/start" "${case_dir}/testresults"
  cp "${app_dir}/q2p1_fac3d" "${case_dir}/q2p1_fac3d"
  cp "${app_dir}/_data/q2p1_param.dat" "${case_dir}/_data/q2p1_param.dat"
  cp "${app_dir}/_data/MG.dat" "${case_dir}/_data/MG.dat"
  cp "${app_dir}/start/data.TXT" "${case_dir}/start/data.TXT"
  cp "${app_dir}/start/sampleRigidBody.xml" \
    "${case_dir}/start/sampleRigidBody.xml"
  ln -s "${app_dir}/_adc" "${case_dir}/_adc"
  ln -s "${partition_root}/${mesh_folder}" "${case_dir}/_mesh/${mesh_folder}"

  sed -i \
    -e "s|^SimPar@MeshFolder.*|SimPar@MeshFolder = \"${mesh_folder}\"|" \
    -e "s|^SimPar@ProjectFile.*|SimPar@ProjectFile = \"${project_file}\"|" \
    -e 's|^SimPar@ProtocolFile.*|SimPar@ProtocolFile = "_data/prot.txt"|' \
    -e 's|^SimPar@MinMeshLevel.*|SimPar@MinMeshLevel = 1|' \
    -e "s|^SimPar@MaxMeshLevel.*|SimPar@MaxMeshLevel = ${mesh_level}|" \
    -e 's|^SimPar@TimeStep.*|SimPar@TimeStep = 0.001d0|' \
    -e "s|^SimPar@MaxNumStep.*|SimPar@MaxNumStep = ${steps}|" \
    -e 's|^SimPar@BackUpFreq.*|SimPar@BackUpFreq = 1000000|' \
    -e 's|^SimPar@OutputFreq.*|SimPar@OutputFreq = 1000000.0d0|' \
    -e 's|^SimPar@OutputLevel.*|SimPar@OutputLevel = 1|' \
    -e 's|^Velo@MGMinLev.*|Velo@MGMinLev = 2|' \
    -e 's|^Velo@MGMedLev.*|Velo@MGMedLev = 2|' \
    -e "s|^Pres@MGCrsSolverType.*|Pres@MGCrsSolverType = ${solver_type}|" \
    "${case_dir}/_data/q2p1_param.dat"

  {
    echo "backend=${backend}"
    echo "ranks=${ranks}"
    echo "workers=$((ranks - 1))"
    echo "solver_type=${solver_type}"
    echo "repetition=${repetition}"
    echo "mesh_level=${mesh_level}"
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
        (
          cd "${case_dir}"
          /usr/bin/time -p mpirun --bind-to core -np "${ranks}" ./q2p1_fac3d
        ) > "${case_dir}/run.log" 2>&1
        grep -E 'HYPRE_TIMING|BenchForce:|Force:|stopped|error|failed' \
          "${case_dir}/run.log" > "${case_dir}/summary.log" || true
        if [[ "${keep_output}" != 1 ]]; then
          find "${case_dir}/_vtk" "${case_dir}/_dump" -type f -delete
        fi
      done
    done
  done
done

python3 "${repo_root}/tools/cluster_scripts/summarize_hypre_gpu_stage1.py" \
  "${result_root}" --warmup-calls "${warmup_calls}" \
  --output "${result_root}/summary.tsv"

echo "Stage 1 FAC3D results: ${result_root}"
