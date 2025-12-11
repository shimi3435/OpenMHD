#!/bin/bash
#------- qsub option -----------
#PBS -P NIFS25KISC015
#PBS -N 2D_KH_PS
#PBS -q B_M
#PBS -l select=4:ncpus=24:mpiprocs=1:ompthreads=24
#PBS -l walltime=5:00
#PBS -j oe

set -euo pipefail

#------- Program execution -----------
module load openmpi/5.0.7/rocm6.3.3

export OMP_NUM_THREADS=24

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

for param_file in "${PARAM_FILES[@]}"; do
  if [[ ! -f "${param_file}" ]]; then
    echo "[WARN] Parameter file ${param_file} not found. Skipping." >&2
    continue
  fi

  base_name=$(basename "${param_file}" .nml)
  timestamp=$(date +%Y%m%d-%H%M%S)
  run_dir="data/${base_name}_${timestamp}"
  mkdir -p "${run_dir}"

  cp "${param_file}" "${run_dir}/params.nml"
  echo "[INFO] Starting run ${base_name} -> ${run_dir}"

  mpirun --display-map --map-by NUMA -x UCX_MAX_RNDV_RAILS=4 ./ap.out "${run_dir}/params.nml" "${run_dir}"
done
