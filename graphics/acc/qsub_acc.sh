#!/bin/bash
################################################################
#
# ACC Analysis Submission Wrapper
# Detects machine (Derecho/Casper) and submits the Python job.
#
# Usage: ./qsub_acc.sh
#
################################################################

root=$PWD
machine=${NCAR_HOST}

if [ -z "$machine" ]; then
    machine=$(hostname)
fi

# --- 1. Machine Specific Setup ---
case $machine in
    *"casper"* )
        account=NMMM0015
        queue=casper
        nodes=1
        ncpus=36
        mem=109GB
        walltime=00:30:00 
        priority="economy"
        ;;
    *"derecho"* )
        account=NMMM0043
        queue=develop
        nodes=1
        ncpus=36
        mem=109GB
        walltime=00:30:00
        priority="economy"
        ;;
    *)
        echo "[!] Error: ${machine} not valid."
        exit 1
        ;;
esac

# --- 2. Create the Job Script ---
# This generates a temporary file to be picked up by qsub
cat <<EOF > temp_acc_job.pbs
#!/bin/bash
#PBS -N ACC_Analysis
#PBS -A ${account}
#PBS -q ${queue}
#PBS -l select=1:ncpus=${ncpus}:mpiprocs=${ncpus}:mem=${mem}
#PBS -l walltime=${walltime}
#PBS -l job_priority=${priority}
#PBS -j oe
#PBS -k eod

# Load Environment
source /etc/profile.d/z00_modules.sh
module load conda/latest
conda activate npl
module list
export PYTHONDONTWRITEBYTECODE=1

# Move to the working directory
cd ${root}

echo "----------------------------------------------------------------"
echo "JOB STARTED: \$(date)"
echo "NODE: \$(hostname)"
echo "PYTHON: \$(which python)"
echo "----------------------------------------------------------------"

# Execute the ACC Analysis script
# The Python script reads its configuration from config_acc.yaml
python run_acc_analysis.py

echo "----------------------------------------------------------------"
echo "JOB ENDED: \$(date)"
exit 0
EOF

# --- 3. Submit and Cleanup ---
chmod +x temp_acc_job.pbs
echo "[*] Submitting ACC job to ${queue} on ${machine}..."
job_id=$(qsub temp_acc_job.pbs)

if [ $? -eq 0 ]; then
    echo "[SUCCESS] Job submitted: ${job_id}"
else
    echo "[ERROR] Submission failed."
fi

exit 0
