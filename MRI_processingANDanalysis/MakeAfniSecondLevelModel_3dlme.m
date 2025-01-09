clear all

%%%Variables to make script%%%
%%%variable names%%%
part_ids = [3:44]; %[3:10,12:15,17:20,23:28]
removeSubs = [3,25,31,36,40]; %remove 5 (2 for missing and 2 for movement,1 for both
%removeSubs = [3,36,40]; %remove 3 (3 for movement)
%removeSubs = [25,31,36]; %remove 3 (3 for missing)
TestID = 4; %1 = onesample ttest, 2 = paired sample ttest, ...
            % 3 = covariate analysis, 4 = LME
IncludeBetweenSubj = 1; %only use if you are including a covariate. ...
                % Check the covariate index variable below to know what your options are
mask_dsat_yn = 0; %0 = no mask, 1 = mask data. If yes, then check ...
                  % mask name and mask dir option below to maek sure they are accurate.
clustCalc_yn = 1; %you will need to have more than 4 particpants to do this....
                  % you can only do this if you are using 3dttest++
LoadBehFile = 1; %if you want to load a behavioral file to use, set this to 1;


%check that all of these are correct
ContrastPrefix = 'Cont_by_Word_lme';
behFile = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/backup_data/FinalDataFile.csv';
data_dir = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/derivatives/afni/';
mask_dir = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/derivatives/afni/SecondLevelModels/FullMod_Blur6_Polort5_3to44_4SubsRem/LME_wordBYcong/';
SLM_Fold_name = 'SecondLevelModels';
FLM_Fold_name = 'FirstLevelModels';
SLM_ModFold_name = 'FullMod_Blur6_Polort5_3to44_5SubsRem';
FLM_ModFold_name = 'FullMod_Blur6_Polort5';
maskName = 'mask_intersection2+tlrc';
TestName = [{'OneSampleTTest'},{'PairedSampleTTest'},{'CovariateAnalysis'},{'LME'}];
slash = '\'; 
DQ = '"';



behData = readtable([behFile],"FileType","text");
currentTime = datetime(now,'ConvertFrom','datenum');
ContrastName = ['LME_wordBYcong'];
prefix = ['CongByWord_interaction'];


behData.DscoreG = behData.Dscore<0;
filename = [ContrastName,'.sh'];
%dirCov =[data_dir,SLM_Fold_name,SLM_ModFold_name,'/',ContrastName,'/covs.txt'];


if ~exist([data_dir,SLM_Fold_name,'/',SLM_ModFold_name,'/',ContrastName], 'dir')
    mkdir([data_dir,SLM_Fold_name,'/',SLM_ModFold_name,'/',ContrastName])
end

cd ([data_dir, SLM_Fold_name,'/',SLM_ModFold_name,'/',ContrastName,'/'])

if IncludeBetweenSubj == 0
    lmeModel = 'Cong*Word*Dscore';qVar = 'Dscore';ranEff = '~1';sstype = '3';
    %DT_colnames = [{'Subj'},{'Cong'},{'Word'},{'Dscore'}];
    DT_colnames = [{'Subj'},{'Cong'},{'Word'},{'Dscore'}];


%%make afni single group ttest script%%%
fileID = fopen(filename, 'w');

%%%check if the file was opened successfully%%%
if fileID == -1
    error('Unable to open the file for writing.')
end

%%%write text and variables to the file%%%
fprintf(fileID, ['module load afni\n\n']);
fprintf(fileID, ['#!/bin/tcsh -xef\n\n']);
fprintf(fileID, '#creation date: %s\n\n',currentTime);
fprintf(fileID, '# ---------------------- set process variables ----------------------\n\n');
fprintf(fileID, 'export dirA=%s%s/%s/\n\n',data_dir,FLM_Fold_name,FLM_ModFold_name);
fprintf(fileID, 'export results_dir=%s%s/%s/%s/\n\n',data_dir,SLM_Fold_name,SLM_ModFold_name,ContrastName);
fprintf(fileID, 'cd $results_dir\n\n');

fprintf(fileID, '# ---------------------- process the data ----------------------\n');
fprintf(fileID, '3dLME -prefix %s%s\n',prefix, slash);
fprintf(fileID, '          -model %s%s%s%s\n',DQ,lmeModel,DQ,slash);
fprintf(fileID, '          -qVars %s%s%s%s\n',DQ,qVar,DQ,slash);
fprintf(fileID, '          -ranEff %s%s%s%s\n',DQ,ranEff,DQ,slash);
fprintf(fileID, '          -SS_type %s%s%s%s\n',DQ,sstype,DQ,slash);
fprintf(fileID, '          -dataTable %s\n',slash);
fprintf(fileID, '%s  %s  %s  %s  InputFile %s\n',DT_colnames{1}, DT_colnames{2},DT_colnames{3},DT_colnames{4},slash);
for i = 1:size(part_ids,2)
%% 
%% 
    i
    if ismember(behData.Subj(i),removeSubs)==0
        fprintf(fileID, '%d    Incongruent    Death    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[4]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Death    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[7]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    Life    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[16]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Life    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[13]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    Me    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[25]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Me    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[22]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    NotMe    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[31]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    NotMe    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[34]" %s\n',...
            behData.Subj(i),behData.Dscore(i),behData.Subj(i),behData.Subj(i),slash);


        % fprintf(fileID, '%d    Incongruent    Death    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[4]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Congruent    Death    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[7]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Incongruent    Life    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[16]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Congruent    Life    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[13]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Incongruent    Me    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[25]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Congruent    Me    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[22]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Incongruent    NotMe    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[31]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);
        % fprintf(fileID, '%d    Congruent    NotMe    %d    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[34]" %s\n',...
        %     behData.Subj(i),behData.DscoreG(i),behData.Subj(i),behData.Subj(i),slash);

    
    else
    end

end
else IncludeBetweenSubj == 1
    lmeModel = 'Cong*Word';ranEff = '~1';sstype = '3';
    %DT_colnames = [{'Subj'},{'Cong'},{'Word'},{'Dscore'}];
    DT_colnames = [{'Subj'},{'Cong'},{'Word'}];
    %%make afni single group ttest script%%%
    fileID = fopen(filename, 'w');

%%%check if the file was opened successfully%%%
if fileID == -1
    error('Unable to open the file for writing.')
end

%%%write text and variables to the file%%%
fprintf(fileID, ['module load afni\n\n']);
fprintf(fileID, ['#!/bin/tcsh -xef\n\n']);
fprintf(fileID, '#creation date: %s\n\n',currentTime);
fprintf(fileID, '# ---------------------- set process variables ----------------------\n\n');
fprintf(fileID, 'export dirA=%s%s/%s/\n\n',data_dir,FLM_Fold_name,FLM_ModFold_name);
fprintf(fileID, 'export results_dir=%s%s/%s/%s/\n\n',data_dir,SLM_Fold_name,SLM_ModFold_name,ContrastName);
fprintf(fileID, 'cd $results_dir\n\n');

fprintf(fileID, '# ---------------------- process the data ----------------------\n');
fprintf(fileID, '3dLME -prefix %s%s\n',prefix, slash);
fprintf(fileID, '          -model %s%s%s%s\n',DQ,lmeModel,DQ,slash);
fprintf(fileID, '          -ranEff %s%s%s%s\n',DQ,ranEff,DQ,slash);
fprintf(fileID, '          -SS_type %s%s%s%s\n',DQ,sstype,DQ,slash);
fprintf(fileID, '          -dataTable %s\n',slash);
fprintf(fileID, '%s  %s  %s  InputFile %s\n',DT_colnames{1}, DT_colnames{2},DT_colnames{3},slash);
for i = 1:size(part_ids,2)
%% 
%% 
    i
    if ismember(behData.Subj(i),removeSubs)==0
        fprintf(fileID, '%d    Incongruent    Death   "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[4]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Death     "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[7]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    Life    "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[16]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Life      "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[13]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    Me      "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[25]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    Me        "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[22]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Incongruent    NotMe   "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[31]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);
        fprintf(fileID, '%d    Congruent    NotMe     "$dirA/sub-%d/stats.sub-%d_run-01+tlrc[34]" %s\n',...
            behData.Subj(i),behData.Subj(i),behData.Subj(i),slash);

    else
    end

end
end

fclose(fileID);


command = ['chmod +rwx ',char(39),filename,char(39)];system(command);


