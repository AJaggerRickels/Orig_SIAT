#!/bin/bash
#SBATCH --job-name=afni${SUBJECT}
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00

#Test that $STUDY is defined
if [ -z "$STUDY" ]; then echo "STUDY is not defined."; exit 1; fi
if [ -z "$SUBJECT" ]; then echo "SUBJECT is not defined."; exit 1; fi
if [ -z "$NUMRUNS" ]; then echo "NUMRUNS is not defined."; exit 1; fi

echo "STUDY   = " $STUDY
echo "SUBJECT =  " $SUBJECT
echo "NUMRUNS =  " $NUMRUNS
echo "BLURSZ = " $BLURSZ
echo "POLORTSZ = " $POLORTSZ
echo "FOLDER = " $FOLDER

module load afni

# Prepare some writeable bind-mount points.

#AFNI_HOST_CACHE=$HOME/.cache/afni
#mkdir -p ${FMRIPREP_HOST_CACHE}

# WORK_DIR="$STUDY/work"
# mkdir -p ${WORK_DIR}

# Prepare derivatives folder
# mkdir -p ${STUDY}/derivatives/afni
 mkdir -p ${STUDY}/derivatives/afni/FirstLevelModels/${FOLDER}/sub-${SUBJECT}/stimuli
 cp ${STUDY}/derivatives/afni/regressors/sub-${SUBJECT}/* ${STUDY}/derivatives/afni/FirstLevelModels/${FOLDER}/sub-${SUBJECT}/stimuli

export SSL_NO_VERIFY=1
export CURL_CA_BUNDLE=""
#export SINGULARITYENV_TEMPLATEFLOW_HOME="/templateflow"
export SINGULARITYENV_SSL_NO_VERIFY=1
export SINGULARITYENV_CURL_CA_BUNDLE=""
# export SINGULARITYENV_FREESURFERHOME=$FREESURFER_HOME
# export LOCAL_FREESURFER_DIR=$FREESURFER_HOME
# mkdir -p ${LOCAL_FREESURFER_DIR}
# Remove IsRunning files from FreeSurfer
# find ${LOCAL_FREESURFER_DIR}/sub-$SUBJECT/ -name "*IsRunning*" -type f -delete
# cp -a $FREESURFER_HOME/license.txt $LOCAL_FREESURFER_DIR/

SINGULARITY_CMD="/usr/bin/singularity"

export FLM_DIR=${STUDY}/derivatives/afni/FirstLevelModels/${FOLDER}/sub-${SUBJECT}
export FMRIPrep_DIR=${STUDY}/derivatives/fmriprep/sub-${SUBJECT}/ses-1/func
export SOURCE_DIR=${STUDY}/sourcedata/sub-${SUBJECT}/ses-1/func
export STIM_DIR=${FLM_DIR}/stimuli

cd ${FLM_DIR}

	if ((${NUMRUNS}==2))
	then
		for ((nr=1;nr<=$NUMRUNS;nr++))
		do
			echo "$nr"
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz ${FLM_DIR}/
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz ${FLM_DIR}/
			
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_desc-confounds_timeseries.tsv ${FLM_DIR}/


			3dmerge -1blur_fwhm ${BLURSZ} -doall -prefix sub-${SUBJECT}-run-0${nr}_blur.nii \
			 sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz


			3dTstat -prefix rm.mean_sub-${SUBJECT}-run-0${nr}.nii sub-${SUBJECT}-run-0${nr}_blur.nii

			3dcalc -a sub-${SUBJECT}-run-0${nr}_blur.nii -b rm.mean_sub-${SUBJECT}-run-0${nr}.nii \
			 -c sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz \
			 -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
			 -prefix sub-${SUBJECT}-run-0${nr}_scale.nii

			3dmask_tool -inputs *mask.nii.gz -union -prefix full_mask.nii
			for reg in trans_x trans_y trans_z rot_x rot_y rot_z; do
			    awk -v col=$reg 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_desc-confounds_timeseries.tsv > sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp.txt;
		    	sed '1d' sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp.txt > sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp_h.txt
		    	
			done
			3dDeconvolve -input sub-${SUBJECT}-run-0${nr}_scale.nii                            \
		    -mask sub-${SUBJECT}_ses-1_task-siat_run-0${nr}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz						     \
		    -polort ${POLORTSZ}                                                                \
		    -num_stimts 20  								\
		    -local_times  								\
		    -GOFORIT 3  								\
			-stim_times 1 stimuli/sub-${SUBJECT}_InstructOnsets_run-0${nr}.txt 'BLOCK(10,1)'      \
			-stim_label 1 Instructions_Onsets                                        \
			-stim_times 2 stimuli/sub-${SUBJECT}_DeathWords_Me_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 2 DeathWords_Me                                           \
			-stim_times 3 stimuli/sub-${SUBJECT}_DeathWords_Other_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 3 DeathWords_Other                                      \
			-stim_times 4 stimuli/sub-${SUBJECT}_DeathWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 4 DeathWords_Prac                                      \
			-stim_times 5 stimuli/sub-${SUBJECT}_LifeWords_Me_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 5 LifeWords_Me                                        \
			-stim_times 6 stimuli/sub-${SUBJECT}_LifeWords_Other_run-0${nr}.txt 'BLOCK(2,1)'       \
			-stim_label 6 LifeWords_Other                                       \
			-stim_times 7 stimuli/sub-${SUBJECT}_LifeWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 7 LifeWords_Prac                                        \
			-stim_times 8 stimuli/sub-${SUBJECT}_MeWords_Life_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 8 MeWords_Life                                \
			-stim_times 9 stimuli/sub-${SUBJECT}_MeWords_Death_run-0${nr}.txt 'BLOCK(2,1)'       \
			-stim_label 9 MeWords_Death                                     \
			-stim_times 10 stimuli/sub-${SUBJECT}_MeWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 10 MeWords_Prac                                     \
			-stim_times 11 stimuli/sub-${SUBJECT}_OthersWords_Life_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 11 OthersWords_Life                                     \
			-stim_times 12 stimuli/sub-${SUBJECT}_OthersWords_Death_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 12 OthersWords_Death                                     \
			-stim_times 13 stimuli/sub-${SUBJECT}_OthersWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 13 OthersWords_Prac                                   \
			-stim_times 14 stimuli/sub-${SUBJECT}_InCorr_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 14 InCorr                                                       \
		    -stim_file 15 sub-${SUBJECT}-trans_x_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 15 -stim_label 15 trans_x_01   \
		    -stim_file 16 sub-${SUBJECT}-trans_y_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 16 -stim_label 16 trans_y_01  \
		    -stim_file 17 sub-${SUBJECT}-trans_z_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 17 -stim_label 17 trans_z_01    \
		    -stim_file 18 sub-${SUBJECT}-rot_x_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 18 -stim_label 18 rot_x_01     \
		    -stim_file 19 sub-${SUBJECT}-rot_y_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 19 -stim_label 19 rot_y_01     \
		    -stim_file 20 sub-${SUBJECT}-rot_z_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 20 -stim_label 20 rot_z_01     \
		    -gltsym 'SYM: +DeathWords_Me -DeathWords_Other'					\
		    -glt_label	1 DeathWords_Me-DeathWords_Other					\
		    -gltsym 'SYM: +LifeWords_Other -LifeWords_Me'			\
		    -glt_label 2 LifeWords_Other-LifeWords_Me					\
		    -gltsym 'SYM: +MeWords_Death -MeWords_Life'					\
		    -glt_label 3 MeWords_Death-MeWords_Life						\
		    -gltsym 'SYM: +OthersWords_Life -OthersWords_Death'				\
	    	    -glt_label	4 OthersWords_Life-OthersWords_Death					\
		    -gltsym 'SYM: +0.25*DeathWords_Me +0.25*LifeWords_Other +0.25*MeWords_Death +0.25*OthersWords_Life -0.25*DeathWords_Other -0.25*LifeWords_Me -0.25*MeWords_Life -0.25*OthersWords_Death'				\
	    	    -glt_label	5 Incong-Cong                                                          \
	            -fout -tout -x1D X_run-0${nr}.xmat.1D -xjpeg X_run-0${nr}.jpg                                  \
		    -x1D_uncensored X_run-0${nr}.nocensor.xmat.1D                                       \
		    -fitts fitts.sub-${SUBJECT}_run-0${nr}                                                       \
		    -errts errts.sub-${SUBJECT}_run-0${nr}                                                     \
		    -bucket stats.sub-${SUBJECT}_run-0${nr}
		    tcsh -x stats.REML_cmd
		done
				for reg in trans_x trans_y trans_z rot_x rot_y rot_z; do
			cat sub-${SUBJECT}-${reg}_run-01_all_tmp_h.txt sub-${SUBJECT}-${reg}_run-02_all_tmp_h.txt > sub-${SUBJECT}-${reg}_allRuns_h.txt
		done
	 	3dDeconvolve -input *_scale.nii                           \
	    -mask full_mask.nii						     \
	    -polort ${POLORTSZ}                                                                \
	    -num_stimts 20  								\
	    -local_times  								\
	    -GOFORIT 3  								\
		-stim_times 1 stimuli/sub-${SUBJECT}_InstructOnsets_allRuns.txt 'BLOCK(10,1)'      \
		-stim_label 1 Instructions_Onsets                                        \
		-stim_times 2 stimuli/sub-${SUBJECT}_DeathWords_Me_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 2 DeathWords_Me                                           \
		-stim_times 3 stimuli/sub-${SUBJECT}_DeathWords_Other_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 3 DeathWords_Other                                      \
		-stim_times 4 stimuli/sub-${SUBJECT}_DeathWords_prac_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 4 DeathWords_Prac                                      \
		-stim_times 5 stimuli/sub-${SUBJECT}_LifeWords_Me_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 5 LifeWords_Me                                        \
		-stim_times 6 stimuli/sub-${SUBJECT}_LifeWords_Other_allRuns.txt 'BLOCK(2,1)'       \
		-stim_label 6 LifeWords_Other                                       \
		-stim_times 7 stimuli/sub-${SUBJECT}_LifeWords_prac_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 7 LifeWords_Prac                                        \
		-stim_times 8 stimuli/sub-${SUBJECT}_MeWords_Life_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 8 MeWords_Life                                \
		-stim_times 9 stimuli/sub-${SUBJECT}_MeWords_Death_allRuns.txt 'BLOCK(2,1)'       \
		-stim_label 9 MeWords_Death                                     \
		-stim_times 10 stimuli/sub-${SUBJECT}_MeWords_prac_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 10 MeWords_Prac                                     \
		-stim_times 11 stimuli/sub-${SUBJECT}_OthersWords_Life_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 11 OthersWords_Life                                     \
		-stim_times 12 stimuli/sub-${SUBJECT}_OthersWords_Death_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 12 OthersWords_Death                                     \
		-stim_times 13 stimuli/sub-${SUBJECT}_OthersWords_prac_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 13 OthersWords_Prac                                   \
		-stim_times 14 stimuli/sub-${SUBJECT}_InCorr_allRuns.txt 'BLOCK(2,1)'      \
		-stim_label 14 InCorr                                                       \
	    -stim_file 15 sub-${SUBJECT}-trans_x_allRuns_h.txt'[0]' -stim_base 15 -stim_label 15 trans_x_01   \
	    -stim_file 16 sub-${SUBJECT}-trans_y_allRuns_h.txt'[0]' -stim_base 16 -stim_label 16 trans_y_01  \
	    -stim_file 17 sub-${SUBJECT}-trans_z_allRuns_h.txt'[0]' -stim_base 17 -stim_label 17 trans_z_01    \
	    -stim_file 18 sub-${SUBJECT}-rot_x_allRuns_h.txt'[0]' -stim_base 18 -stim_label 18 rot_x_01     \
	    -stim_file 19 sub-${SUBJECT}-rot_y_allRuns_h.txt'[0]' -stim_base 19 -stim_label 19 rot_y_01     \
	    -stim_file 20 sub-${SUBJECT}-rot_z_allRuns_h.txt'[0]' -stim_base 20 -stim_label 20 rot_z_01     \
	    -gltsym 'SYM: +DeathWords_Me -DeathWords_Other'					\
	    -glt_label	1 DeathWords_Me-DeathWords_Other					\
	    -gltsym 'SYM: +LifeWords_Other -LifeWords_Me'			\
	    -glt_label 2 LifeWords_Other-LifeWords_Me					\
	    -gltsym 'SYM: +MeWords_Death -MeWords_Life'					\
	    -glt_label 3 MeWords_Death-MeWords_Life						\
	    -gltsym 'SYM: +OthersWords_Life -OthersWords_Death'				\
		-glt_label	4 OthersWords_Life-OthersWords_Death					\
	    -gltsym 'SYM: +0.25*DeathWords_Me +0.25*LifeWords_Other +0.25*MeWords_Death +0.25*OthersWords_Life -0.25*DeathWords_Other -0.25*LifeWords_Me -0.25*MeWords_Life -0.25*OthersWords_Death'\	
	    -glt_label	5 Incong-Cong								\
	    -fout -tout -x1D X_allRuns.xmat.1D -xjpeg X_allRuns.jpg                                  \
	    -x1D_uncensored X_allRuns.nocensor.xmat_allRuns.1D                                       \
	    -fitts fitts.sub-${SUBJECT}_allRuns                                                      \
	    -errts errts.sub-${SUBJECT}_allRuns                                                     \
	    -bucket stats.sub-${SUBJECT}_allRuns 

	   tcsh -x stats.REML_cmd
else
			nr=1
			echo "$nr"
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz ${FLM_DIR}/
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz ${FLM_DIR}/
			
			cp ${FMRIPrep_DIR}/sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_desc-confounds_timeseries.tsv ${FLM_DIR}/


			3dmerge -1blur_fwhm ${BLURSZ} -doall -prefix sub-${SUBJECT}-run-0${nr}_blur.nii \
			 sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz


			3dTstat -prefix rm.mean_sub-${SUBJECT}-run-0${nr}.nii sub-${SUBJECT}-run-0${nr}_blur.nii

			3dcalc -a sub-${SUBJECT}-run-0${nr}_blur.nii -b rm.mean_sub-${SUBJECT}-run-0${nr}.nii \
			 -c sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz \
			 -expr 'c * min(200, a/b*100)*step(a)*step(b)' \
			 -prefix sub-${SUBJECT}-run-0${nr}_scale.nii

			3dmask_tool -inputs *mask.nii.gz -union -prefix full_mask.nii
			for reg in trans_x trans_y trans_z rot_x rot_y rot_z; do
			    awk -v col=$reg 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}} print $c} NR>1{print $c}' sub-${SUBJECT}_ses-1_task-siat_run-0"$nr"_desc-confounds_timeseries.tsv > sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp.txt;
		    	sed '1d' sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp.txt > sub-${SUBJECT}-${reg}_run-0${nr}_all_tmp_h.txt
		    	
			done
			3dDeconvolve -input sub-${SUBJECT}-run-0${nr}_scale.nii                            \
		    -mask sub-${SUBJECT}_ses-1_task-siat_run-0${nr}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz						     \
		    -polort ${POLORTSZ}                                                                \
		    -num_stimts 20  								\
		    -local_times  								\
		    -GOFORIT 3  								\
			-stim_times 1 stimuli/sub-${SUBJECT}_InstructOnsets_run-0${nr}.txt 'BLOCK(10,1)'      \
			-stim_label 1 Instructions_Onsets                                        \
			-stim_times 2 stimuli/sub-${SUBJECT}_DeathWords_Me_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 2 DeathWords_Me                                           \
			-stim_times 3 stimuli/sub-${SUBJECT}_DeathWords_Other_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 3 DeathWords_Other                                      \
			-stim_times 4 stimuli/sub-${SUBJECT}_DeathWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 4 DeathWords_Prac                                      \
			-stim_times 5 stimuli/sub-${SUBJECT}_LifeWords_Me_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 5 LifeWords_Me                                        \
			-stim_times 6 stimuli/sub-${SUBJECT}_LifeWords_Other_run-0${nr}.txt 'BLOCK(2,1)'       \
			-stim_label 6 LifeWords_Other                                       \
			-stim_times 7 stimuli/sub-${SUBJECT}_LifeWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 7 LifeWords_Prac                                        \
			-stim_times 8 stimuli/sub-${SUBJECT}_MeWords_Life_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 8 MeWords_Life                                \
			-stim_times 9 stimuli/sub-${SUBJECT}_MeWords_Death_run-0${nr}.txt 'BLOCK(2,1)'       \
			-stim_label 9 MeWords_Death                                     \
			-stim_times 10 stimuli/sub-${SUBJECT}_MeWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 10 MeWords_Prac                                     \
			-stim_times 11 stimuli/sub-${SUBJECT}_OthersWords_Life_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 11 OthersWords_Life                                     \
			-stim_times 12 stimuli/sub-${SUBJECT}_OthersWords_Death_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 12 OthersWords_Death                                     \
			-stim_times 13 stimuli/sub-${SUBJECT}_OthersWords_prac_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 13 OthersWords_Prac                                   \
			-stim_times 14 stimuli/sub-${SUBJECT}_InCorr_run-0${nr}.txt 'BLOCK(2,1)'      \
			-stim_label 14 InCorr                                                       \
		    -stim_file 15 sub-${SUBJECT}-trans_x_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 15 -stim_label 15 trans_x_01   \
		    -stim_file 16 sub-${SUBJECT}-trans_y_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 16 -stim_label 16 trans_y_01  \
		    -stim_file 17 sub-${SUBJECT}-trans_z_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 17 -stim_label 17 trans_z_01    \
		    -stim_file 18 sub-${SUBJECT}-rot_x_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 18 -stim_label 18 rot_x_01     \
		    -stim_file 19 sub-${SUBJECT}-rot_y_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 19 -stim_label 19 rot_y_01     \
		    -stim_file 20 sub-${SUBJECT}-rot_z_run-0${nr}_all_tmp_h.txt'[0]' -stim_base 20 -stim_label 20 rot_z_01     \
		    -gltsym 'SYM: +DeathWords_Me -DeathWords_Other'					\
		    -glt_label	1 DeathWords_Me-DeathWords_Other					\
		    -gltsym 'SYM: +LifeWords_Other -LifeWords_Me'			\
		    -glt_label 2 LifeWords_Other-LifeWords_Me					\
		    -gltsym 'SYM: +MeWords_Death -MeWords_Life'					\
		    -glt_label 3 MeWords_Death-MeWords_Life						\
		    -gltsym 'SYM: +OthersWords_Life -OthersWords_Death'				\
	    	    -glt_label	4 OthersWords_Life-OthersWords_Death					\
		    -gltsym 'SYM: +0.25*DeathWords_Me +0.25*LifeWords_Other +0.25*MeWords_Death +0.25*OthersWords_Life -0.25*DeathWords_Other -0.25*LifeWords_Me -0.25*MeWords_Life -0.25*OthersWords_Death'				\
	    	    -glt_label	5 Incong-Cong                                                          \
	            -fout -tout -x1D X_run-0${nr}.xmat.1D -xjpeg X_run-0${nr}.jpg                                  \
		    -x1D_uncensored X_run-0${nr}.nocensor.xmat.1D                                       \
		    -fitts fitts.sub-${SUBJECT}_run-0${nr}                                                       \
		    -errts errts.sub-${SUBJECT}_run-0${nr}                                                     \
		    -bucket stats.sub-${SUBJECT}_run-0${nr}
		    tcsh -x stats.REML_cmd
		
fi

	


	

