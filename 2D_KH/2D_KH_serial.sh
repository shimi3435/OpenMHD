#!/bin/bash
#------- qsub option -----------
#PBS -N 2D_KH_serial
#PBS -j oe

set -euo pipefail

#------- Program execution -----------

cd "${PBS_O_WORKDIR}"

# A space separated list of parameter files can be supplied through PARAM_LIST.
# Otherwise, we look for params_*.nml under params_cases/, and fall back to params.nml.
declare -a PARAM_FILES
if [[ -n "${PARAM_LIST:-}" ]]; then
  read -r -a PARAM_FILES <<< "${PARAM_LIST}"
elif compgen -G "params_cases/*.nml" > /dev/null; then
  mapfile -t PARAM_FILES < <(ls params_cases/*.nml)
else
  PARAM_FILES=("params.nml")
fi

if [[ ${#PARAM_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No parameter files found." >&2
  exit 1
fi

if [[ -n "${RUN_DIR:-}" && ${#PARAM_FILES[@]} -ne 1 ]]; then
  echo "[ERROR] RUN_DIR requires a single parameter file." >&2
  exit 1
fi

for param_file in "${PARAM_FILES[@]}"; do
  if [[ ! -f "${param_file}" ]]; then
    echo "[WARN] Parameter file ${param_file} not found. Skipping." >&2
    continue
  fi

  base_name=$(basename "${param_file}" .nml)
  timestamp=$(date +%Y%m%d_%H%M%S)
  if [[ -n "${RUN_DIR:-}" ]]; then
    run_dir="${RUN_DIR}"
    if [[ "${run_dir}" != /* ]]; then
      run_dir="${PBS_O_WORKDIR}/${run_dir}"
    fi
  else
    run_dir="data/${base_name}_${timestamp}"
  fi
  mkdir -p "${run_dir}"

  cp "${param_file}" "${run_dir}/params.nml"
  echo "[INFO] Starting run ${base_name} -> ${run_dir}"

  ./a.out "${run_dir}/params.nml" "${run_dir}"
done
