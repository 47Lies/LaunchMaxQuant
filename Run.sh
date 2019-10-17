#PBS -l mem=10gb
#PBS -l walltime=1800:00:00
#PBS -j oe

source /bioinfo/users/${LOGNAME}/.bashrc
WD=/data/users/ltaing/DATA_TMP/ltaing/IsoAndSpe
if [ "$HOSTNAME" == "u900-bdd-1-185n-6925.curie.fr" ]
then
 LocalProfile="standard"
else
 LocalProfile="cluster"
fi

cd ${WD}
nextflow run Main.nf -c Nextflow.Config.yaml -profile ${LocalProfile} -resume


