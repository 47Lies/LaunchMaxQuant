#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=240:00:00
#PBS -j oe

source ~/.bashrc
conda activate mono

Mass_Spec_SampleDescriptionFile=SampleDescription.txt
DATABASE=uniprot_sprot.9606.fasta


if [[ ! -z ${PBS_JOBID} ]];
then
cd ${PBS_O_WORKDIR}
WORKDIR=/local/scratch/${PBS_JOBID}
mkdir -p ${WORKDIR}/RAW
#cp RAW/* ${WORKDIR}/RAW
#cp ${DATABASE} ${WORKDIR}/${DATABASE}
WORKERS=${PBS_NUM_PPN}
else
WORKDIR=${PBS_O_WORKDIR}
WORKERS=8
fi
echo "workdir $WORKDIR"
echo "path $PWD"
${PWD}/TemplateWithFiles.py -c\
 -f ${PWD}/${DATABASE}\
 -t ${PWD}\
 -w ${PWD}\
 -p ${WORKERS}\
 -s ${Mass_Spec_SampleDescriptionFile}\
 -r ${PWD}/RAW\
 -m mqpar.xml\
 -o Temporary.mqpar.xml
mono MaxQuant_1.6.14/bin/MaxQuantCmd.exe Temporary.mqpar.xml
mkdir -p RESULTS
mv ${PWD}/fixedCombinedFolder/combined/txt/* RESULTS/.
