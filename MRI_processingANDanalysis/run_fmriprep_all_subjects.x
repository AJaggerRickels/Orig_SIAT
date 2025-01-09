#!/bin/bash
#To run on the cluster do:
#
#ssh nerve1    (if not already logged in)
#cd to working directory

#Usage ./run_fmriprep_all_subjects.x study_dir

if [ "$1" == "" ]; then
	echo "Usage: ./run_fmriprep_all_subjects.x study_dir" 
	exit 1
fi

export STUDY=$1

cd $STUDY

#Process participants.tsv file
#Skip first row, Print first column, remove sub-, remove leading 0, remove blank lines at end.
#cat participants.tsv | awk 'NR > 1 { print$1 }' | sed 's/sub-//' | sed 's/^0*//'
cat participants.tsv | grep "sub-" | awk '{ print $1 }' | sed 's/sub-//' | sed 's/^0*//' > tmp_participants

participants=$(cat tmp_participants)
for participant in $participants
do
	echo "Submitting participant = $participant"
	export SUBJECT=$participant
	sbatch /cm/shared/apps/fmriprep/run_fmriprep_single_subject.sbatch
done

rm -f tmp_participants

