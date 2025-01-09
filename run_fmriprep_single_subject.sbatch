#!/bin/bash
#SBATCH --job-name=fmriprep$SUBJECT
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00

#Test that $STUDY is defined
if [ -z "$STUDY" ]; then echo "STUDY is not defined."; exit 1; fi
if [ -z "$SUBJECT" ]; then echo "SUBJECT is not defined."; exit 1; fi

echo "STUDY   = " $STUDY
echo "SUBJECT =  " $SUBJECT

module load freesurfer

# Prepare some writeable bind-mount points.
TEMPLATEFLOW_HOST_HOME=$HOME/.cache/templateflow
FMRIPREP_HOST_CACHE=$HOME/.cache/fmriprep
mkdir -p ${TEMPLATEFLOW_HOST_HOME}
mkdir -p ${FMRIPREP_HOST_CACHE}

WORK_DIR="$STUDY/work"
mkdir -p ${WORK_DIR}

# Prepare derivatives folder
mkdir -p ${STUDY}/derivatives/fmriprep-22.1.0 

export SSL_NO_VERIFY=1
export CURL_CA_BUNDLE=""
export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
export SINGULARITYENV_SSL_NO_VERIFY=1
export SINGULARITYENV_CURL_CA_BUNDLE=""
export SINGULARITYENV_FREESURFERHOME=$FREESURFER_HOME
export LOCAL_FREESURFER_DIR=$FREESURFER_HOME
mkdir -p ${LOCAL_FREESURFER_DIR}
# Remove IsRunning files from FreeSurfer
find ${LOCAL_FREESURFER_DIR}/sub-$SUBJECT/ -name "*IsRunning*" -type f -delete
cp -a $FREESURFER_HOME/license.txt $LOCAL_FREESURFER_DIR/

SINGULARITY_CMD="/usr/bin/singularity"

$SINGULARITY_CMD run --cleanenv -B $STUDY/sourcedata:/sourcedata -B $STUDY/derivatives/fmriprep-22.1.0:/derivatives -B ~/.cache/templateflow:/templateflow -B $STUDY/work:/work -B $LOCAL_FREESURFER_DIR /cm/shared/apps/fmriprep/fmriprep-22.1.0.simg /sourcedata /derivatives participant --participant-label $SUBJECT -w /work/ -vv --omp-nthreads 8 --nthreads 16 --mem_mb 30000 --output-spaces MNI152NLin2009cAsym:res-2 anat fsnative fsaverage5 --use-aroma --fs-subjects-dir  $STUDY/derivatives/fmriprep-22.1.0/freesurfer-7.2/ --fs-license-file $LOCAL_FREESURFER_DIR/license.txt

