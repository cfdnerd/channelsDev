#!/bin/bash
#SBATCH --job-name=clusterKeepAlive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=61
#SBATCH --ntasks=61
#SBATCH --output=%x-%j.Out
#SBATCH --partition=hpc

set -euo pipefail

echo "[$(date '+%Y-%m-%d %H:%M:%S')] keep-alive job started on ${SLURM_JOB_NODELIST:-unknown}"
echo "This job intentionally performs minimal work to keep a 61-task allocation active."

# Keep the allocation alive with negligible load across all task slots.
echo "Launching ${SLURM_NTASKS:-61} keep-alive processes..."
exec srun --ntasks="${SLURM_NTASKS:-61}" --ntasks-per-node=61 sleep infinity
