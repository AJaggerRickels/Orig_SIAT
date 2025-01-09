#!/bin/tcsh
#script written for use with APE AAT to generate the average blur estimate.

set wkdir = /quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/derivatives/afni/FirstLevelModels/FullMod_Blur6_Polort5


#for testing.
#set subjects = (s21183)
#echo subjects

#set subjects = (3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44)

set subjects = (4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 26 27 28 29 30 32 33 34 35 37 38 39 41 42 43 44)

set sub_string = "sub-"

#foreach subject ($subj_loop)
foreach subj ($subjects)

set final_subj = "${sub_string}${subj}"

echo "who's that? ...Oh shit, it's ${final_subj}!"

set subj_dir = ${wkdir}/${final_subj}
set acf_dir = ${wkdir}

#get rid of files that have already been made (to avoid duplicate rows within files)
cd $subj_dir

if ( -f ${final_subj}_run-01.3dFWHMx.ACF ) then
    echo output dir "${final_subj}_run-01.3dFWHMx.ACF" already exists
    rm  ${final_subj}_run-01.3dFWHMx.ACF
endif

if ( -f TEMP ) then
    echo output dir "TEMP" already exists
    rm  TEMP
endif

if ( -f ${final_subj}_run-01.blur.errts.1d.ACF) then
    echo output dir "${final_subj}_run-01.blur.errts.1d.ACF" already exists
    rm  ${final_subj}_run-01.blur.errts.1d.ACF
endif


#Now compute the acf parameters for each subject

3dFWHMx -detrend -acf ${final_subj}_run-01.3dFWHMx.ACF 	\
#-mask ${subj_dir}/full_mask+tlrc 	\
${subj_dir}/errts.${final_subj}_run-01+tlrc > TEMP.ACF
tail -n 1 TEMP.ACF >> ${final_subj}_run-01.blur.errts.1D.ACF

cat ${final_subj}_run-01.blur.errts.1D.ACF >> $acf_dir/temp.blurs.1D.ACF

end

cd $acf_dir

3dTstat -mean -prefix all_blurs.1D.ACF_Final_CombinedModel temp.blurs.1D.ACF\'

