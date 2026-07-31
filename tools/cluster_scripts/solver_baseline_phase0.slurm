#!/bin/bash
#SBATCH --job-name=ff-solver-p0
#SBATCH --partition=short
#SBATCH --nodelist=tardis
# The GPU is not used; tardis only accepts jobs that request its gres, and
# pinning tardis keeps the baseline on the same hardware as the Stage 1 runs.
#SBATCH --constraint=ampere
#SBATCH --gres=gpu:ampere:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=384G
#SBATCH --time=02:00:00
#SBATCH --output=stage1-work/slurm-%x-%j.out

# Phase 0 CPU solver baseline (see docs/md_docs/solver_baseline_phase0.md).
#
# Runs the cylinder (q2p1_fc_ext) and undeformed FAC3D (q2p1_fac3d)
# benchmarks with the existing pressure coarse-solver types and records
# FF_TIMING / HYPRE_TIMING per-call solver timings. SSE is intentionally
# not part of this harness. Requires a build configured with
# -DENABLE_SOLVER_TIMING=ON (and -DENABLE_HYPRE_TIMING=ON for types 7/8).

set -euo pipefail

module purge
module load gcc/13.2.0
module load openmpi/options/interface/ethernet openmpi/options/cuda/no openmpi/4.1.6

# backend=gpu runs CUDA-hypre builds: loads the CUDA module and serves the
# workers through CUDA MPS (as in the Stage 1 scripts).
backend="${FF_P0_BACKEND:-cpu}"
if [[ "${backend}" == gpu ]]; then
  module load cuda/12.4
fi

repo_root="${FF_P0_ROOT:-${SLURM_SUBMIT_DIR}}"
build_dir="${FF_P0_BUILD:-${repo_root}/stage1-work/build-matched-cpu}"
apps="${FF_P0_APPS:-cylinder fac3d}"
solver_types="${FF_P0_SOLVER_TYPES:-1 2 3 4 7 8}"
velo_types="${FF_P0_VELO_TYPES:-1}"
velo_minlev="${FF_P0_VELO_MINLEV:-}"
repeats="${FF_P0_REPEATS:-3}"
warmup_calls="${FF_P0_WARMUP_CALLS:-2}"
steps_cylinder="${FF_P0_STEPS_CYLINDER:-10}"
steps_fac3d="${FF_P0_STEPS_FAC3D:-4}"
ranks_cylinder="${FF_P0_RANKS_CYLINDER:-4 8}"
ranks_fac3d="${FF_P0_RANKS_FAC3D:-4}"
fac3d_level="${FF_P0_FAC3D_LEVEL:-4}"
keep_output="${FF_P0_KEEP_OUTPUT:-0}"
run_id="${FF_P0_RUN_ID:-phase0-${SLURM_JOB_ID}}"
result_root="${repo_root}/stage1-work/results/${run_id}"

app_binary() {
  case "$1" in
    cylinder) echo "q2p1_fc_ext" ;;
    fac3d) echo "q2p1_fac3d" ;;
    *) echo "Unknown app: $1" >&2; return 1 ;;
  esac
}

app_project() {
  case "$1" in
    cylinder) echo "_adc/2D_FAC/2Dbench.prj" ;;
    fac3d) echo "_adc/3D_FAC_FBM/file.prj" ;;
  esac
}

app_mesh_prefix() {
  case "$1" in
    cylinder) echo "NEWFAC" ;;
    fac3d) echo "FAC3D" ;;
  esac
}

app_steps() {
  case "$1" in
    cylinder) echo "${steps_cylinder}" ;;
    fac3d) echo "${steps_fac3d}" ;;
  esac
}

app_ranks() {
  case "$1" in
    cylinder) echo "${ranks_cylinder}" ;;
    fac3d) echo "${ranks_fac3d}" ;;
  esac
}

for app in ${apps}; do
  binary="$(app_binary "${app}")"
  app_dir="${build_dir}/applications/${binary}"
  if [[ ! -x "${app_dir}/${binary}" ]]; then
    echo "Missing executable: ${app_dir}/${binary}" >&2
    exit 2
  fi
  if [[ ! -f "${app_dir}/$(app_project "${app}")" ]]; then
    echo "Missing project $(app_project "${app}") under ${app_dir}" >&2
    exit 2
  fi
done

mkdir -p "${result_root}/runs" "${result_root}/metadata"

{
  echo "run_id=${run_id}"
  echo "host=$(hostname)"
  echo "job_id=${SLURM_JOB_ID}"
  echo "build_dir=${build_dir}"
  echo "apps=${apps}"
  echo "solver_types=${solver_types}"
  echo "velo_types=${velo_types}"
  echo "velo_minlev=${velo_minlev}"
  echo "repeats=${repeats}"
  echo "warmup_calls=${warmup_calls}"
  echo "steps_cylinder=${steps_cylinder}"
  echo "steps_fac3d=${steps_fac3d}"
  echo "ranks_cylinder=${ranks_cylinder}"
  echo "ranks_fac3d=${ranks_fac3d}"
  echo "fac3d_level=${fac3d_level}"
  echo "backend=${backend}"
  echo "git_commit=$(git -C "${repo_root}" rev-parse HEAD)"
  echo "mpi=$(mpirun --version | head -1)"
} > "${result_root}/metadata/run.env"

if [[ "${backend}" == gpu ]]; then
  export CUDA_VISIBLE_DEVICES=0
  nvidia-smi -q > "${result_root}/metadata/nvidia-smi.txt"
  mps_pipe="${TMPDIR:-/tmp}/ff-p0-mps-${SLURM_JOB_ID}"
  mps_log="${TMPDIR:-/tmp}/ff-p0-mps-log-${SLURM_JOB_ID}"
  mkdir -p "${mps_pipe}" "${mps_log}"
  export CUDA_MPS_PIPE_DIRECTORY="${mps_pipe}"
  export CUDA_MPS_LOG_DIRECTORY="${mps_log}"
  nvidia-cuda-mps-control -d
  nvidia-smi dmon -s pucvmet -d 1 \
    > "${result_root}/metadata/nvidia-smi-dmon.log" 2>&1 &
  dmon_pid=$!
  cleanup() {
    kill "${dmon_pid}" 2>/dev/null || true
    echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true
  }
  trap cleanup EXIT
fi

for app in ${apps}; do
  binary="$(app_binary "${app}")"
  app_dir="${build_dir}/applications/${binary}"
  project_file="$(app_project "${app}")"
  mesh_prefix="$(app_mesh_prefix "${app}")"
  # PyPartitioner dlopens libmetis.so from the application directory; not
  # every application target stages it.
  if [[ ! -f "${app_dir}/libmetis.so" ]]; then
    cp "${build_dir}/extern/libraries/metis-5.1.0/libmetis/libmetis.so" \
      "${app_dir}/libmetis.so"
  fi
  for ranks in $(app_ranks "${app}"); do
    workers=$((ranks - 1))
    mesh_folder="${mesh_prefix}_np${ranks}"
    mesh_file="${app_dir}/_mesh/${mesh_folder}/sub0001/GRID0001.tri"
    if [[ ! -f "${mesh_file}" ]]; then
      (
        cd "${app_dir}"
        LD_LIBRARY_PATH="${app_dir}:${LD_LIBRARY_PATH:-}" \
          python3 "${repo_root}/tools/PyPartitioner.py" \
          "${workers}" 1 1 "${mesh_folder}" "${project_file}"
      ) > "${result_root}/metadata/partition-${app}-np${ranks}.log" 2>&1 || true
    fi
    if [[ ! -f "${mesh_file}" ]]; then
      echo "Partitioning failed for ${app} np${ranks}; see" \
        "${result_root}/metadata/partition-${app}-np${ranks}.log" >&2
      exit 3
    fi
  done
done

stage_case() {
  local app=$1
  local ranks=$2
  local solver_type=$3
  local velo_type=$4
  local repetition=$5
  local case_dir=$6
  local binary app_dir project_file mesh_folder steps

  binary="$(app_binary "${app}")"
  app_dir="${build_dir}/applications/${binary}"
  project_file="$(app_project "${app}")"
  mesh_folder="$(app_mesh_prefix "${app}")_np${ranks}"
  steps="$(app_steps "${app}")"

  mkdir -p "${case_dir}/_data" "${case_dir}/_dump" "${case_dir}/_gmv" \
    "${case_dir}/_mesh" "${case_dir}/_ns" "${case_dir}/_vtk" \
    "${case_dir}/solution" "${case_dir}/start" "${case_dir}/testresults"
  cp "${app_dir}/${binary}" "${case_dir}/${binary}"
  cp "${app_dir}/_data/q2p1_param.dat" "${case_dir}/_data/q2p1_param.dat"
  cp "${app_dir}/_data/MG.dat" "${case_dir}/_data/MG.dat"
  cp "${app_dir}/start/data.TXT" "${case_dir}/start/data.TXT"
  cp "${app_dir}/start/sampleRigidBody.xml" \
    "${case_dir}/start/sampleRigidBody.xml"
  ln -s "${app_dir}/_adc" "${case_dir}/_adc"
  ln -s "${app_dir}/_mesh/${mesh_folder}" "${case_dir}/_mesh/${mesh_folder}"

  sed -i \
    -e "s|^SimPar@MeshFolder.*|SimPar@MeshFolder = \"${mesh_folder}\"|" \
    -e 's|^SimPar@ProtocolFile.*|SimPar@ProtocolFile = "_data/prot.txt"|' \
    -e "s|^SimPar@MaxNumStep.*|SimPar@MaxNumStep = ${steps}|" \
    -e 's|^SimPar@BackUpFreq.*|SimPar@BackUpFreq = 1000000|' \
    -e 's|^SimPar@OutputFreq.*|SimPar@OutputFreq = 1000000.0d0|' \
    -e "s|^Pres@MGCrsSolverType.*|Pres@MGCrsSolverType = ${solver_type}|" \
    -e "s|^Velo@MGCrsSolverType.*|Velo@MGCrsSolverType = ${velo_type}|" \
    "${case_dir}/_data/q2p1_param.dat"

  if [[ "${app}" == fac3d ]]; then
    sed -i \
      -e "s|^SimPar@ProjectFile.*|SimPar@ProjectFile = \"${project_file}\"|" \
      -e 's|^SimPar@MinMeshLevel.*|SimPar@MinMeshLevel = 1|' \
      -e "s|^SimPar@MaxMeshLevel.*|SimPar@MaxMeshLevel = ${fac3d_level}|" \
      -e 's|^SimPar@TimeStep.*|SimPar@TimeStep = 0.001d0|' \
      -e 's|^SimPar@ApplyFAC3DMeshDeformation.*|SimPar@ApplyFAC3DMeshDeformation = No|' \
      -e 's|^SimPar@OutputLevel.*|SimPar@OutputLevel = 1|' \
      -e 's|^Velo@MGMinLev.*|Velo@MGMinLev = 2|' \
      -e 's|^Velo@MGMedLev.*|Velo@MGMedLev = 2|' \
      "${case_dir}/_data/q2p1_param.dat"
  fi

  # Optional override of the velocity coarse level (e.g. to align it with
  # Pres@MGMedLev, required by the velocity MUMPS solver).
  if [[ -n "${velo_minlev}" ]]; then
    sed -i \
      -e "s|^Velo@MGMinLev.*|Velo@MGMinLev = ${velo_minlev}|" \
      -e "s|^Velo@MGMedLev.*|Velo@MGMedLev = ${velo_minlev}|" \
      "${case_dir}/_data/q2p1_param.dat"
  fi

  {
    echo "app=${app}"
    echo "backend=${backend}"
    echo "ranks=${ranks}"
    echo "workers=$((ranks - 1))"
    echo "solver_type=${solver_type}"
    echo "velo_type=${velo_type}"
    echo "repetition=${repetition}"
    echo "steps=${steps}"
  } > "${case_dir}/case.env"
}

for app in ${apps}; do
  binary="$(app_binary "${app}")"
  for ranks in $(app_ranks "${app}"); do
    for solver_type in ${solver_types}; do
     for velo_type in ${velo_types}; do
      for repetition in $(seq 1 "${repeats}"); do
        name="${app}-${backend}-np${ranks}-type${solver_type}-vt${velo_type}-rep${repetition}"
        case_dir="${result_root}/runs/${name}"
        stage_case "${app}" "${ranks}" "${solver_type}" "${velo_type}" \
          "${repetition}" "${case_dir}"
        run_exit_code=0
        (
          cd "${case_dir}"
          /usr/bin/time -p mpirun --bind-to core -np "${ranks}" "./${binary}"
        ) > "${case_dir}/run.log" 2>&1 || run_exit_code=$?
        echo "run_exit_code=${run_exit_code}" >> "${case_dir}/case.env"
        grep -E 'FF_TIMING|HYPRE_TIMING|BenchForce:|Force:|stopped|error|failed' \
          "${case_dir}/run.log" > "${case_dir}/summary.log" || true
        if [[ "${keep_output}" != 1 ]]; then
          find "${case_dir}/_vtk" "${case_dir}/_dump" -type f -delete
        fi
      done
     done
    done
  done
done

python3 "${repo_root}/tools/cluster_scripts/summarize_hypre_gpu_stage1.py" \
  "${result_root}" --warmup-calls "${warmup_calls}" \
  --output "${result_root}/summary.tsv"

echo "Phase 0 baseline results: ${result_root}"
