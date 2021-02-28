#PBS -l nodes=1:ppn=60
#PBS -l mem=200gb
#PBS -l walltime=240:00:00
#PBS -q 'mpi'
#PBS -j oe

#source ~/.bashrc
#conda activate mono

Mass_Spec_SampleDescriptionFile=SampleDescription.txt





if [[ ! -z ${PBS_JOBID} ]];
then
cd ${PBS_O_WORKDIR}
WORKDIR=/local/scratch/${PBS_JOBID}
mkdir -p ${WORKDIR}/RAW
cp RAW/*.raw ${WORKDIR}/RAW
cp ${DATABASE} ${WORKDIR}/${DATABASE}
WORKERS=${PBS_NUM_PPN}
else
WORKDIR=`pwd -P`
WORKERS=1
fi
echo $WORKDIR
${WORKDIR}/TemplateWithFiles.py -c\
 -f ${WORKDIR}/${DATABASE}\
 -t ${WORKDIR}\
 -w ${WORKDIR}\
 -p ${WORKERS}\
 -s ${Mass_Spec_SampleDescriptionFile}\
 -r ${WORKDIR}/RAW\
 -m mqpar.xml\
 -o Temporary.mqpar.xml
mono MaxQuant_1.6.14/bin/MaxQuantCmd.exe Temporary.mqpar.xml --dryrun
mkdir -p RESULTS
mv ${WORKDIR}/fixedCombinedFolder/combined/txt/* RESULTS/.
