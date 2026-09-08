#!/usr/bin/env bash

# Build and install FeatFloWer, prepare a gendie or heat case, and print the
# command needed to run it from any directory.
#
# Examples:
#   ./cmake_cheat_setup_from_Raphael.sh -M gendie /scratch/cases/die_02
#   ./cmake_cheat_setup_from_Raphael.sh -M heat -i /scratch/4508/CASE /scratch/4508/RUN
#
# Optional environment variables:
#   INSTALL_PREFIX=/path/to/install  Installation directory (default: ../INSTALL)
#   BUILD_DIR=/path/to/build         Build directory (default: build-<module>)
#   BUILD_JOBS=4                     Parallel Ninja jobs
#   MPI_RANKS=4                      MPI ranks shown in the run command
#   INPUT_SOURCE=/path/to/input      Same as -i/--input

set -Eeuo pipefail

usage() {
  cat >&2 <<EOF
Usage: $0 [-M gendie|heat] [-i INPUT] <runtime-case>

Options:
  -M, --module NAME  Required simulation module: gendie or heat
  -i, --input DIR    Input template copied into the runtime case
  -h, --help         Show this help

With -M gendie, INPUT defaults to the installed RH_LOC example and is copied
to <runtime-case>/e3d_input. With -M heat, INPUT defaults to the installed HEAT
test case and is copied to <runtime-case>/heat_input. A custom heat INPUT must
contain heat.s3d, sampleRigidBody.xml, top-level OFF geometry, and optionally
meshDir. An already populated heat_input is preserved when INPUT is omitted.
EOF
}

simulation_module=
input_source=${INPUT_SOURCE:-}
case_arg=

while [[ $# -gt 0 ]]; do
  case "$1" in
    -M|--module)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; usage; exit 2; }
      simulation_module=$2
      shift 2
      ;;
    --module=*)
      simulation_module=${1#*=}
      shift
      ;;
    -i|--input)
      [[ $# -ge 2 ]] || { echo "Error: $1 requires a value." >&2; usage; exit 2; }
      input_source=$2
      shift 2
      ;;
    --input=*)
      input_source=${1#*=}
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      if [[ $# -ne 1 || -n "${case_arg}" ]]; then
        echo "Error: specify exactly one runtime case folder." >&2
        usage
        exit 2
      fi
      case_arg=$1
      shift
      ;;
    -*)
      echo "Error: unknown option '$1'." >&2
      usage
      exit 2
      ;;
    *)
      if [[ -n "${case_arg}" ]]; then
        echo "Error: specify exactly one runtime case folder." >&2
        usage
        exit 2
      fi
      case_arg=$1
      shift
      ;;
  esac
done

if [[ $# -gt 0 || -z "${case_arg}" ]]; then
  usage
  exit 2
fi
if [[ -z "${simulation_module}" ]]; then
  echo "Error: select a simulation module with -M gendie or -M heat." >&2
  usage
  exit 2
fi
if [[ "${simulation_module}" != gendie && "${simulation_module}" != heat ]]; then
  echo "Error: module must be 'gendie' or 'heat', not '${simulation_module}'." >&2
  exit 2
fi

for command_name in cmake ninja python3 realpath; do
  if ! command -v "${command_name}" >/dev/null 2>&1; then
    echo "Error: required command '${command_name}' was not found." >&2
    exit 1
  fi
done

case_dir=$(realpath -m -- "${case_arg}")

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir=${BUILD_DIR:-"${script_dir}/build-${simulation_module}"}
install_prefix=${INSTALL_PREFIX:-"${script_dir}/../INSTALL"}
build_jobs=${BUILD_JOBS:-4}
mpi_ranks=${MPI_RANKS:-4}

if [[ -n "${input_source}" ]]; then
  input_source=$(realpath -m -- "${input_source}")
  if [[ ! -d "${input_source}" ]]; then
    echo "Error: input template directory does not exist: ${input_source}" >&2
    exit 2
  fi
fi

echo "Repository:      ${script_dir}"
echo "Build directory: ${build_dir}"
echo "Install prefix:  ${install_prefix}"
echo "Module:          ${simulation_module}"
echo "Runtime case:    ${case_dir}"
if [[ -n "${input_source}" ]]; then
  echo "Input template:  ${input_source}"
fi

# 1. Configure. Heat and gendie both require the CGAL-enabled application tree.
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

# 4. Prepare the selected module's input and runtime folders.
mkdir -p "${case_dir}"

if [[ "${simulation_module}" == gendie ]]; then
  gendie_dir="${install_prefix}/bin/q2p1_gendie"
  default_input="${gendie_dir}/_ianus/DIE/RH_LOC"
  source_input=${input_source:-"${default_input}"}
  target_input="${case_dir}/e3d_input"
  parameter_template="${gendie_dir}/_data_BU/q2p1_paramV_DIE_test.dat"
  launcher="${gendie_dir}/e3d_start_yaml.py"

  for required in "${source_input}" "${parameter_template}" "${launcher}"; do
    if [[ ! -e "${required}" ]]; then
      echo "Error: required gendie input was not found: ${required}" >&2
      exit 1
    fi
  done

  mkdir -p "${target_input}" "${case_dir}/_data_BU"
  if [[ "$(realpath -m -- "${source_input}")" != "$(realpath -m -- "${target_input}")" ]]; then
    cp -a "${source_input}/." "${target_input}/"
  fi

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

  echo
  echo "Gendie build, installation, and case preparation completed."
  echo "Run the simulation from any directory with:"
  printf 'python3 %q -C %q -f %q -y %q -n %q -x\n' \
    "${launcher}" "${case_dir}" "${target_input}" \
    "${case_dir}/plan_test.yaml" "${mpi_ranks}"
else
  heat_dir="${install_prefix}/bin/heat"
  target_input="${case_dir}/heat_input"
  launcher="${heat_dir}/heat_start.py"
  default_input="${heat_dir}/_ianus/HEAT"

  if [[ -n "${input_source}" ]]; then
    source_input=${input_source}
  elif [[ -f "${target_input}/heat.s3d" && \
          -f "${target_input}/sampleRigidBody.xml" ]]; then
    source_input=${target_input}
  else
    source_input=${default_input}
  fi

  for required in "${launcher}" "${source_input}/heat.s3d" \
                  "${source_input}/sampleRigidBody.xml"; do
    if [[ ! -e "${required}" ]]; then
      echo "Error: required heat input was not found: ${required}" >&2
      echo "Use -i INPUT or install the bundled _ianus/HEAT test case." >&2
      exit 1
    fi
  done

  if [[ "$(realpath -m -- "${source_input}")" != \
        "$(realpath -m -- "${target_input}")" ]]; then
    mkdir -p "${target_input}"
    cp -a "${source_input}/." "${target_input}/"
  fi

  for required in "${launcher}" "${target_input}/heat.s3d" \
                  "${target_input}/sampleRigidBody.xml"; do
    if [[ ! -e "${required}" ]]; then
      echo "Error: required heat input was not found: ${required}" >&2
      echo "Use -i INPUT or populate ${target_input}." >&2
      exit 1
    fi
  done
  if ! find "${target_input}" -maxdepth 1 -type f \
      \( -name '*.off' -o -name '*.OFF' \) -print -quit | grep -q .; then
    echo "Error: no top-level OFF geometry found in ${target_input}." >&2
    exit 1
  fi

  echo
  echo "Heat build, installation, and case preparation completed."
  echo "Run the simulation from any directory with:"
  printf 'python3 %q -C %q -f %q -n %q\n' \
    "${launcher}" "${case_dir}" "${target_input}" "${mpi_ranks}"
fi
