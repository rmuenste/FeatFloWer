#!/usr/bin/env bash

# Build and install FeatFloWer, prepare the RH_LOC DIE example, and print the
# command needed to run it.
#
# Usage:
#   ./cmake_cheat_setup_from_Raphael.sh /scratch/cases/die_02
#
# Optional environment variables:
#   INSTALL_PREFIX=/path/to/install  Installation directory (default: ../INSTALL)
#   BUILD_JOBS=4                     Parallel Ninja jobs
#   MPI_RANKS=4                      MPI ranks shown in the run command

set -Eeuo pipefail

usage() {
  echo "Usage: $0 /scratch/<case-folder>" >&2
}

if [[ $# -ne 1 ]]; then
  usage
  exit 2
fi

for command_name in cmake ninja python3 realpath; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "Error: required command '${command_name}' was not found." >&2
    exit 1
  fi
done

case_dir=$(realpath -m -- "$1")
if [[ "${case_dir}" != /scratch/* ]]; then
  echo "Error: the case folder must be located below /scratch." >&2
  echo "Resolved path: ${case_dir}" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="${script_dir}/build-gendie"
install_prefix=${INSTALL_PREFIX:-"${script_dir}/../INSTALL"}
build_jobs=${BUILD_JOBS:-4}
mpi_ranks=${MPI_RANKS:-4}

echo "Repository:      ${script_dir}"
echo "Build directory: ${build_dir}"
echo "Install prefix:  ${install_prefix}"
echo "Case directory:  ${case_dir}"

# 1. Configure.
cmake -S "${script_dir}" -B "${build_dir}" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_APPLICATIONS=ON \
  -DUSE_CGAL=ON \
  -DUSE_HYPRE=OFF \
  -DUSE_PE=OFF \
  -DENABLE_FBM_ACCELERATION=OFF

# 2. Build the whole configured tree. The install step installs every active
# install rule and therefore requires all corresponding artifacts to exist.
ninja -C "${build_dir}" -j"${build_jobs}"

# 3. Install.
cmake --install "${build_dir}" --prefix "${install_prefix}"

# 4. Prepare the case folder.
gendie_dir="${install_prefix}/bin/q2p1_gendie"
example_input="${gendie_dir}/_ianus/DIE/RH_LOC"
parameter_template="${gendie_dir}/_data_BU/q2p1_paramV_DIE_test.dat"
launcher="${gendie_dir}/e3d_start_yaml.py"

if [[ ! -d "${example_input}" ]]; then
  echo "Error: installed example input was not found: ${example_input}" >&2
  exit 1
fi
if [[ ! -f "${parameter_template}" ]]; then
  echo "Error: installed parameter template was not found: ${parameter_template}" >&2
  exit 1
fi
if [[ ! -f "${launcher}" ]]; then
  echo "Error: installed launcher was not found: ${launcher}" >&2
  exit 1
fi

mkdir -p "${case_dir}/e3d_input" "${case_dir}/_data_BU"
cp -a "${example_input}/." "${case_dir}/e3d_input/"

cat >"${case_dir}/plan_test.yaml" <<'EOF'
options:
  die-simulation: true
  ResolutionLevel: 2

stages:
  init:
    steps:
      - solver: MomentumEquation
        param_file: _data_BU/q2p1_paramV_DIE_test.dat
  final:
    steps:
      - solver: MomentumEquation
        param_file: _data_BU/q2p1_paramV_DIE_test.dat
EOF

# This case-local override is needed because Hypre is disabled above.
sed 's/^Pres@MGCrsSolverType *= *8/Pres@MGCrsSolverType = 1/' \
  "${parameter_template}" \
  >"${case_dir}/_data_BU/q2p1_paramV_DIE_test.dat"

# 5. Print the simulation command instead of starting it. All paths are
# absolute, and -C makes the simulation use the requested case directory, so
# the printed command can be executed from anywhere.
echo
echo "Build, installation, and case preparation completed."
echo "Run the simulation from any directory with:"
printf 'python3 %q -C %q -f %q -y %q -n %q -x\n' \
  "${launcher}" \
  "${case_dir}" \
  "${case_dir}/e3d_input" \
  "${case_dir}/plan_test.yaml" \
  "${mpi_ranks}"
