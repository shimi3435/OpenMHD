#!/bin/bash
# Submit a chain of qsub jobs with varying parameters.
# Usage: ./run_param_sweep.sh [template_nml] [run_list] [job_script]
set -euo pipefail

TEMPLATE=${1:-params.nml}
RUN_LIST=${2:-param_runs.list}
JOB_SCRIPT=${3:-2D_KH_mpi.sh}
PARAM_DIR="params_cases"

if [[ ! -f "${TEMPLATE}" ]]; then
  echo "[ERROR] Template ${TEMPLATE} not found." >&2
  exit 1
fi

if [[ ! -f "${RUN_LIST}" ]]; then
  echo "[ERROR] Run definition file ${RUN_LIST} not found." >&2
  exit 1
fi

if [[ ! -f "${JOB_SCRIPT}" ]]; then
  echo "[ERROR] Job script ${JOB_SCRIPT} not found." >&2
  exit 1
fi

mkdir -p "${PARAM_DIR}"

update_param() {
  local file=$1
  local key=$2
  local value=$3
  python3 - "${file}" "${key}" "${value}" <<'PY'
import sys,re
path,key,value = sys.argv[1:4]
with open(path) as f:
    lines = f.readlines()
pat = re.compile(r'^(\s*' + re.escape(key) + r'\s*=\s*)(.*)', re.IGNORECASE)
for idx,line in enumerate(lines):
    m = pat.match(line)
    if m:
        lines[idx] = f"{m.group(1)}{value}\n"
        break
else:
    raise SystemExit(f"key '{key}' not found in {path}")
with open(path,'w') as f:
    f.writelines(lines)
PY
}

wait_for_job() {
  local job_id=$1
  local interval=${WAIT_INTERVAL:-30}
  echo "[INFO] Waiting for job ${job_id} to finish..."
  while true; do
    if qstat "${job_id}" >/dev/null 2>&1; then
      sleep "${interval}"
    else
      break
    fi
  done
}

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" || "${line#\#}" != "$line" ]] && continue
  IFS=' ' read -r -a parts <<< "$line"
  run_name=${parts[0]}
  if [[ -z "${run_name}" ]]; then
    continue
  fi
  params_file="${PARAM_DIR}/${run_name}.nml"
  cp "${TEMPLATE}" "${params_file}"

  for kv in "${parts[@]:1}"; do
    key=${kv%%=*}
    value=${kv#*=}
    update_param "${params_file}" "${key}" "${value}"
  done

  env_param="PARAM_LIST=${params_file}"
  job_output=$(qsub -v "${env_param}" "${JOB_SCRIPT}")
  job_id=$(echo "${job_output}" | awk '{print $1}')
  echo "[INFO] Submitted ${run_name} as ${job_id}"
  wait_for_job "${job_id}"
done < "${RUN_LIST}"
