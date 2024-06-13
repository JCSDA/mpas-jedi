#!/bin/bash
################################################################
#
## Driving script to submit PBS jobs to conduct combinations of
#  verifications according to the binning method and domain. The
#  Otional:
#  - hydro: True/False
#
# Usage: ./qsub.sh
#
################################################################

root=$PWD

machine=`hostname`

case $machine in
    *"casper"* )
        account=NMMM0015
        queue=casper
        nodes=1
        ncpus=36
        mpiprocs=36
        mem=109GB
        ompthreads=2
        walltime=01:00:00
        job_priority=economy
        ;;
    *"derecho"* )
        account=NMMM0043
        queue=main
        nodes=1
        ncpus=128
        mpiprocs=128
        mem=235GB
        ompthreads=1
        walltime=00:45:00
        job_priority=economy
        ;;

    *)
        echo "${machine} not valid"
        ;;
esac
####################### VARIABLES SETUP ########################
rootDir=/glade/derecho/scratch/ivette/pandac/
ctpPath=/glade/derecho/scratch/bjung/saca_new_mesh/obs_processing/convert.mpas_10kmL2/output_fixed_clearCTP/
staticPath=/glade/campaign/mmm/parc/ivette/pandac/codeBuild/MPAS-Model/init_files/

experiments="ColdStart_clean2.1"
experiments+=",CloudDirectInsertion_clean2.1"
experiments+=",3dhybrid-allsky_O60-3kmI60km_clean2.1"

names="ColdStart"
names+=",DIfrom0hr"
names+=",CyclingDA"

binningMethod=("CTPLayering" "CloudFraction")
domains=("fullDisk" "finer" "coarser")
dateIni=2018041500
dateEnd=2018043000
delta=24
fcstHr=18
hydro=true

for binning in ${binningMethod[@]}; do
  for domain in ${domains[@]}; do

  outputFolder=${binning}_${domain}
  if ${hydro}; then
    outputFolder=${binning}_${domain}_hydro
  fi
  mkdir -p ${outputFolder}
  cd ${outputFolder}

cat <<EOF > ${binning}_${domain}.sh

#!/bin/sh
##############################################################
# Script to submit ${binning}_${domain}.sh script
##############################################################
#!/bin/bash

#PBS -N ${binning}_${domain}
#PBS -l select=${nodes}:ncpus=${ncpus}:mpiprocs=${mpiprocs}:mem=${mem}:ompthreads=${ompthreads}
#PBS -l walltime=${walltime}
#PBS -q ${queue}
#PBS -A ${account}
#PBS -l job_priority=${job_priority}
#PBS -j oe
#PBS -o job.out
#PBS -e job.err

#
# load modules:
# ===========
source /etc/profile.d/z00_modules.sh
module load conda/latest
conda activate npl
module list
export PYTHONDONTWRITEBYTECODE=1

cd $PWD

echo "`date` STARTED DATE"

### Run the python script
echo "Binning by ${binning} and ${domain}"
python ${root}/conditionalVerification.py -dir ${rootDir} -ctp ${ctpPath} -static ${staticPath} -exp "${experiments}" -n "${names}" -di ${dateIni} -de ${dateEnd} -dt ${delta} -fhr ${fcstHr} -b ${binning} -hy ${hydro} -r ${domain}

echo "`date` ENDED DATE"
exit 0

EOF

chmod +x ${binning}_${domain}.sh

echo "Submiting job ${binning}_${domain}.sh ..."
qsub ${binning}_${domain}.sh
cd ..

  done
done
exit 0
