#!/bin/bash
# Torque Configuration
#PBS -l walltime=500:00:00
#PBS -l mem=128gb
#PBS -l nodes=1:ppn=32
#PBS -N MxQunt_BgDt
#PBS -j oe


#
#We need to have the . 'dot' as decimal separator. In our current structure the Environment variable is set to french where the decimal separator is ','
#
#
# The script take a sample description file
#  Matrix file
#   -FileName
#   -Experiment
#   -Fraction
#   -PTMS
#   -ParamGroupIndices
#  *One line per FileName
#  *Typical description of MaxQuant sample option
#
#
# 2 Steps
#  Step I) Build a mqpar.xml that will be used by in house python script
#	Python 3 /data/users/ltaing/miniconda3/bin/python
#	Temporary file directory in the local scratch -t /local/scratch/${PBS_JOBID}
#		Unable to predict it
#	Fasta file to be searched -f /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/SP-Human-UP000005640-012018.fasta
#	The working directory -w /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/
#	The number of threads to be used -p 32
#	The sample description file -s /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/SampleDescription.txt
#	The Termofisher mass spectromery .raw files directory; All .Raw files in a single directory -r /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/RAW
#	The maxquant parameter xml template -m /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/mqpar.template.1.6.2.6.xml
#	The output parameter file -o /data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/mqpar.Silac.AllSample.xml
#  Step II) Run MaxQuant with the newly created mqpar.xml
#
#

export SINGULARITYENV_LANG="en_US.UTF-8"
export LANG="en_US.UTF-8"

########################################################################
#
#
#                        *** Coniguration Part ***
#
#
########################################################################

########################################################################
#	Where Everything is
#	What you smust change
#
WD=/bioinfo/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5

########################################################################
#
#	Your Files location
#	What you should change
#
#	Sample description of your files
#	Database fasta file of the tested proteins
#	Thermo Fischer Raw files directory
#
#
SampleDescriptionFile=${WD}/SampleDescription.txt
DatabaseFastaFile=${WD}/SP-Human-UP000005640-012018.fasta
RAW_Files_DIRECTORY=${WD}/RAW

cd ${WD}
########################################################################
#
#	Executables and Parameters files location
#	What you should not change unless you known what you are doing
#
#	Template of 
#	Database fasta file of the tested proteins
#	Thermo Fischer Raw files directory
#
#
TemplateMqparFile=/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/mqpar.template.1.6.2.6.xml
MxQuntCmd=/bioinfo/users/ltaing/SOFTS/MaxQuantVersion1.6.2.6/bin/MaxQuantCmd.exe
MonoSIMG=/bioinfo/users/ltaing/mono-5.16.0.220.simg
TemporaryMqparFile=${WD}/mqpar.Silac.AllSample.xml


if [[ -f "$TemplateMqparFile" ]];
then
 singularity exec ${MonoSIMG} mono ${MxQuntCmd} --create ${TemplateMqparFile}
fi

/data/users/ltaing/miniconda3/bin/python scripts/TemplateWithFiles.py -c\
 -t /local/scratch/${PBS_JOBID}\
 -f ${DatabaseFastaFile}\
 -w ${WD}\
 -p 32\
 -s ${SampleDescriptionFile}\
 -r ${RAW_Files_DIRECTORY}\
 -m ${TemplateMqparFile}\
 -o ${TemporaryMqparFile}
singularity exec ${MonoSIMG} mono ${MxQuntCmd} ${TemporaryMqparFile}
