#!/bin/bash
#To run on the cluster do:
#
#ssh nerve1    (if not already logged in)
#cd to working directory

#Usage ./run_fmriprep_all_subjects.x study_dir blursize polortsize foldername

if [ "$1" == "" ]; then
	echo "Usage: ./run_fmriprep_all_subjects.x study_dir" 
	exit 1
fi
if [ "$2" == "" ]; then
	echo "Usage: ./run_fmriprep_all_subjects.x study_dir" 
	exit 1
fi
if [ "$3" == "" ]; then
	echo "Usage: ./run_fmriprep_all_subjects.x study_dir" 
	exit 1
fi
if [ "$4" == "" ]; then
	echo "Usage: ./run_fmriprep_all_subjects.x study_dir" 
	exit 1
fi

export STUDY=$1
export BLURSZ=$2
export POLORTSZ=$3
export FOLDER=$4

cd $STUDY

#Process participants.tsv file
#Skip first row, Print first column, remove sub-, remove leading 0, remove blank lines at end.
#cat participants.tsv | awk 'NR > 1 { print$1 }' | sed 's/sub-//' | sed 's/^0*//'
cat participants.tsv | grep "sub-" | awk '{ print $1 }' | sed 's/sub-//' | sed 's/^0*//' > tmp_participants
participants=$(cat tmp_participants)

export count=0
for participant in $participants
do
	export count=`expr $count + 1`
	cat participants.tsv |  sed '1d' | sed -n "${count}p" |awk '{ print $3 }'  > tmp_siat_runs
	echo "Submitting participant = $participant"
	export SUBJECT=$participant
	export NUMRUNS=$(cat tmp_siat_runs)
	echo ${NUMRUNS}
	echo "Output data to " $STUDY "/derivatives/afni/FirstLevelModels/" $FOLDER
	
	sbatch /quobyte/nerve-data/dsalat/2389_trt_main/arickels/Code/FirstLevelAnalysis/run_afni_single_subject_FullModel_newcont.sbatch
done

rm -f tmp_participants
rm -f tmp_siat_runs

