clear all
%%This regressor maker is for the SIAT updated on 06/17/2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%set up parameters and file info for the script. Check to make sure%%%%%%
%%%and make sure that all of this info matches with your current file%%%%%%
%%%structure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subjects = [3:44];
% if you have not done fMRIprep yet, you can only run a section of this 
% script. Select 0 if you HAVE NOT run fMRIprep and 1 if you HAVE run fMRIprep
fMRIprepDone = 1; 
%AFter you run fMRIpre, this will make a copy of the events file and move
%it to the derivitves folder
    cp_source2deriv =1;
    %after you run fMRIprep, this will make afni stimuli as well as give you
    %the option to clacuate the d-score and reaction time csv files.
    MakeAFNIstim = 1;
        %if makeAFNIstim is set to 1, you can either only make afni regressors
        %(MakeRegressors == 1) or you can calcuate the dscore and reaction times by
        %setting MakeRegressors to 0 and CalcDscore and CalcRTS to 1. 
        MakeRegressors =0;
        CalcDScore = 1;
        CalcRTS = 1;
%set up the dscore and RTS file names
DS_filename = 'Dscore_All42Subs'
RTS_filename = 'RTS_All42Subs'
%you dont need to chagne these, leave them as is
fmriprepFolder = ['fmriprep/sub-'];
afniFolder = 'afni/';

%Right now the script is only set up to run the "full model" so leave this
%as is for now,
RegressorType = "Full_Model"
%set these directories as follows
%This is the location of the raw data both behavioral and MRI
RawDataFolder = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/backup_data/BehavData/';
%This is the location of your BIDs organzied data (you have already run the
%DCM2BIDS scripts
SourcedataFolder = ['/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/sourcedata/sub-'];
%This is the location of your preprocessed fMRIprep data
DerivativesFolder = ['/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/derivatives/'];
%This is the location that you would like the calcuated dscores and RTS to
%be set too
BehavDataFolder = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/backup_data/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Check from this point UP that all this information is correct before
      %running this script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DScores = [];
RTS = [];

if fMRIprepDone ==0; 
  %this block will create an events.tsv file for the siat
  for S = 1:size(Subjects,2)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section sets up variable names
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subID = num2str(Subjects(S)); %subject id
    cd (RawDataFolder)
    fileNames = dir (['sub-',subID,'-free*']);
    blockNames = [{'Death|Life'},{'Me|NotMe'},{'Prac_Death/Me|Life/NotMe'},{'Death/Me|Life/NotMe'},...
        {'Life|Death'},{'Prac_Life/Me|Death/NotMe'},{'Life/Me|Death/NotMe'};{'Life|Death'},...
        {'Me|NotMe'},{'Prac_Life/Me|Death/NotMe'},{'Life/Me|Death/NotMe'},{'Death|Life'},...
        {'Prac_Death/Me|Life/NotMe'},{'Death/Me|Life/NotMe'}];
    blockWords = [{'Die'},{'Funeral'},{'Lifeless'},{'Suicide'},{'Deceased'};...
        {'Live'},{'Thrive'},{'Breathing'},{'Alive'},{'Survive'};...
        {'Myself'},{'Mine'},{'I'},{'My'},{'Self'};...
        {'Them'},{'Other'},{'They'},{'Their'},{'Theirs'}];
    
    VariableNames = [{'onset'},{'version'},{'duration'},{'run'},{'trial#'},...
            {'responsetime'},{'response'},...
            {'block'},{'stimtype'},{'stim_id'},{'corrvsincor'},{'blockname'}, {'stim'}];
    
    RespOpts = [-1,1;1,1;-1,2;1,2;-1,3;1,3;-1,4;1,4];
    CorRespInx.Vers1 = {{[1,4]};{[5,8]};{[1,4,5,8]};{[1,4,5,8]};{[2,3]};{[2,3,5,8]};{[2,3,5,8]}};
    CorRespInx.Vers2 = {{[2,3]};{[5,8]};{[2,3,5,8]};{[2,3,5,8]};{[1,4]};{[1,4,5,8]};{[1,4,5,8]}};
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section does a quick check about the number of runs and the
      %version of the SIat that was completed.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for r = 1:size(fileNames,1)
        runNum = ['run',num2str(r)];
        fileInfo.(runNum)= detectImportOptions(fileNames(r,1).name, "Range","A1");
        t.(runNum)= readtable(fileNames(r,1).name,"Range","A2"); %first run of SIAT
        t.(runNum)(:,8:end) = [];
        timeOut(r,1) = {fileInfo.(runNum).VariableOptions(1,6).Name(2:end)};
        if fileInfo.(runNum).VariableOptions(1,8).Name((end-12):(end-9)) == 'run1'
            dispOut(r,1) = {[runNum,' was Version 1']};
        elseif fileInfo.(runNum).VariableOptions(1,8).Name((end-12):(end-9)) == 'run2'
            dispOut(r,1) = {[runNum,' was Version 2']};
        else
        end
    end
    
     if size(fileNames,1) ==1
         disp(['There was only one run and it was: ',dispOut{1,1}(end-9:end)]);
     elseif size(fileNames,1) ==2
         if dispOut{1,1}(end) ==dispOut{2,1}(end)
             disp(['There were two runs and both were the same version: ',dispOut{1,1}(end-9:end)]);
         else
             disp(['There were two runs and they were different versions: ',dispOut{1,1}(end-9:end),' then ',dispOut{2,1}(end-9:end)]);
         end
         if str2num(timeOut{1,1}(end-3:end))<str2num(timeOut{2,1}(end-3:end))
             disp(['Run1 was ran before Run2, files were named correctly'])
         else str2num(timeOut{1,1}(end-3:end))>str2num(timeOut{2,1}(end-3:end))
             disp(['Run2 was ran before Run1, files were named incorrectly, PLEASE CHECK YOUR FILE NAMES!'])
         end
     else
     end
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section creates a matrix/table for EACH RUN of the SIAT that
      %includes variables of interests from the raw data file. 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     for r = 1:size(fileNames,1) %for each run of the SIAT
        
        runNum = ['run',num2str(r)];
        WordOnset_indx(r,1) = {find(t.(runNum).event_code == 5)};
        WordRespons_indx(r,1)= {find(t.(runNum).event_code == 6)};
        Instructions_indx(r,1) = {find(t.(runNum).event_code == 7)};
        x = [t.(runNum).absolute_time(WordOnset_indx{r,1}),(ones(size(WordOnset_indx{r,1},1),1)+1),ones(size(WordOnset_indx{r,1},1),1),...
            t.(runNum).trial_number(WordOnset_indx{r,1}),zeros(size(WordOnset_indx{r,1},1),1),...
            zeros(size(WordOnset_indx{r,1},1),1)];
        y = [ones(size(WordRespons_indx{r,1},1),1),t.(runNum).trial_number(WordRespons_indx{r,1}),str2double(t.(runNum).response_time(WordRespons_indx{r,1})),str2double(t.(runNum).response(WordRespons_indx{r,1}))];
        x_str = [cellstr(num2str(ones(size(WordOnset_indx{r,1},1),1))),t.(runNum).trial_type(WordOnset_indx{r,1}),t.(runNum).result(WordOnset_indx{r,1})];
        Instruct = [ones(size(Instructions_indx{r,1},1),1),t.(runNum).absolute_time(Instructions_indx{r,1})];
    
    clear WordOnset WordRespons_indx Instructions_indx    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section cleans up nans
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:length(x) 
            if isempty(find(x(i,4)==y(y(:,1)==1,2)));
                x(i,5:6) = NaN;
            else
                x(i,5:6) = y(find(x(i,4)==y(y(:,1)==1,2)),3:4);
            end
        end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section sets up the variables for each of the 7 blocks of the
      %siat, and adds a column to code for the verison of the siat
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:7
            blockIndx = find(contains(x_str(:,2),num2str(j)));
            x(blockIndx,7) = j;
            if dispOut{1,1}(end-8:end) == 'Version 1';
               x_str(blockIndx,2) = blockNames(1,j);
               for w = 1:size(blockIndx)
                    isWord = cellfun(@(x)strcmp(x,x_str(blockIndx(w),3)),blockWords);
                    [row,col] = find(isWord);
                    x(blockIndx(w),8) = row;
                    x(blockIndx(w),9) = col;
                    m = find(RespOpts(:,1) == x(blockIndx(w),6) &RespOpts(:,2)==row);
                    m2 = find(CorRespInx.Vers1{j}{1}==m);
                    if isempty(m2)==1;
                        m2 = 0;
                    else; 
                        m2 =1;
                    end;
                    x(blockIndx(w),10) = m2;
                    
               end
               
            else dispOut{1,1}(end-8:end) == 'Version 2';
               x_str(blockIndx,2) = blockNames(2,j);
               for w = 1:size(blockIndx)
                    isWord = cellfun(@(x)strcmp(x,x_str(blockIndx(w),3)),blockWords);
                    [row,col] = find(isWord);
                    x(blockIndx(w),8) = row;
                    x(blockIndx(w),9) = col;
                    m = find(RespOpts(:,1) == x(blockIndx(w),6) &RespOpts(:,2)==row);
                    m2 = find(CorRespInx.Vers2{j}{1}==m);
                    if isempty(m2)==1;
                        m2 = 0;
                    else; 
                        m2 =1;
                    end;
                    x(blockIndx(w),10) = m2;
               
               end
               
            end
         
        end
        if dispOut{1,1}(end-8:end) == 'Version 1'; 
            x = [x(:,1),ones(length(x),1),x(:,2:end)];
        else dispOut{1,1}(end-8:end) == 'Version 2';
            x = [x(:,1),(ones(length(x),1)+1),x(:,2:end)];
        end
        Runs(r,1) = {table(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7),...
            x(:,8),x(:,9),x(:,10),x(:,11),x_str(:,2),x_str(:,3))};
        Runs{r,1}.Properties.VariableNames = VariableNames;
        Inst(r,1) = {Instruct(:,2)};
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %this section sets up some directories, creates some additional
      %variables, and saves the tsv (in bids format for fMRIprep), 
      %and saves csv, and txt files for future anlaysis with AFNI
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    if ~exist([RawDataFolder,'working/sub-',subID,'/'], 'dir');    mkdir([RawDataFolder,'working/sub-',subID,'/']);    end
    if ~exist([SourcedataFolder,subID,'/ses-1/func/'], 'dir');    mkdir([SourcedataFolder,subID,'/ses-1/func/']);    end
   
    cd([SourcedataFolder,subID,'/ses-1/func/'])
        x((Runs{r,1}.stimtype==1),(end+1)) =1 ;
        x((Runs{r,1}.stimtype==2),(end+1)) =1 ;
        x((Runs{r,1}.stimtype==3),(end+1)) =1 ;
        x((Runs{r,1}.stimtype==4),(end+1)) =1 ;
        x(:,(end+1)) =(Runs{r,1}.corrvsincor.*Runs{r,1}.stimtype==1) ;
        x(:,(end+1)) =(Runs{r,1}.corrvsincor.*Runs{r,1}.stimtype==2) ;
        x(:,(end+1)) =(Runs{r,1}.corrvsincor.*Runs{r,1}.stimtype==3) ;
        x(:,(end+1)) =(Runs{r,1}.corrvsincor.*Runs{r,1}.stimtype==4) ;
        Runs{r,1} = [Runs{r,1},table(x(:,(end-7)),x(:,(end-6)),x(:,(end-5)),x(:,(end-4)),x(:,(end-3)),x(:,(end-2)),x(:,(end-1)),x(:,end))]; Runs{r,1}.Properties.VariableNames = [VariableNames,{'deathwords'},{'lifewords'},{'mewords'},{'otherswords'},{'deathwordsCO'},{'lifewordsCO'},{'mewordsCO'},{'otherswordsCO'}];
        %writetable(Runs{r,1},['sub-',subID,'_StimFile_',runNum,'_2.csv']);
        
        stimFile = Runs{r,1};
       
        for i = 1: size(stimFile,1)
            %practice blocks
            if stimFile.block(i) == 1 | stimFile.block(i) == 2 | stimFile.block(i) == 5
                %DeathWords
                if stimFile.deathwords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i) ={['C_Prac_DW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) ==0
                        stimFile.stimLabel(i) ={['IC_Prac_DW_',char(stimFile.stim(i))]};
                    end
                %LifeWords
                elseif stimFile.lifewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Prac_LW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) ==0
                        stimFile.stimLabel(i)={['IC_Prac_LW_',char(stimFile.stim(i))]};
                    end
                %MeWords
                elseif stimFile.mewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Prac_MW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Prac_MW_',char(stimFile.stim(i))]};
                    end
                %NotMeWords
                else stimFile.otherswords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Prac_NW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Prac_NW_',char(stimFile.stim(i))]};
                    end
                end
            elseif stimFile.blockname(i) =="Death/Me|Life/NotMe" |stimFile.blockname(i) =="Prac_Death/Me|Life/NotMe" 
                if stimFile.deathwords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_InCong_DW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_InCong_DW_',char(stimFile.stim(i))]};
                    end
                %LifeWords
                elseif stimFile.lifewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_InCong_LW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_InCong_LW_',char(stimFile.stim(i))]};
                    end
                %MeWords
                elseif stimFile.mewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_InCong_MW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_InCong_MW_',char(stimFile.stim(i))]};
                    end
                %NotMeWords
                else stimFile.otherswords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_InCong_NW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_InCong_NW_',char(stimFile.stim(i))]};
                    end
                end
            else stimFile.blockname(i) =="Life/Me|Death/NotMe" |stimFile.blockname(i) =="Prac_Life/Me|Death/NotMe" 
                if stimFile.deathwords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Cong_DW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Cong_DW_',char(stimFile.stim(i))]};
                    end
                %LifeWords
                elseif stimFile.lifewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Cong_LW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Cong_LW_',char(stimFile.stim(i))]};
                    end
                %MeWords
                elseif stimFile.mewords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Cong_MW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Cong_MW_',char(stimFile.stim(i))]};
                    end
                %NotMeWords
                else stimFile.otherswords(i) ==1
                    if stimFile.corrvsincor(i) == 1
                        stimFile.stimLabel(i)={['C_Cong_NW_',char(stimFile.stim(i))]};
                    else stimFile.corrvsincor(i) == 0
                        stimFile.stimLabel(i)={['IC_Cong_NW_',char(stimFile.stim(i))]};
                    end
                end
            end
        end
        
        Runs{r,1} = stimFile;
        writetable(Runs{r,1},['sub-',subID,'_ses-1_task-siat_run-0',char(string(r)),'_events.tsv'],"FileType","text", 'Delimiter', 'tab', 'WriteVariableNames',true);
        writetable(Runs{r,1},[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_ses-1_task-siat_run-0',char(string(r)),'_events.csv'],"FileType","text", 'Delimiter', ' ');
        writetable(Runs{r,1},['sub-',subID,'_ses-1_task-siat_run-0',char(string(r)),'_events.csv'],"FileType","text", 'Delimiter', ' ');
 
        if dispOut{1,1}(end-8:end) == 'Version 1';
           writematrix(Inst{r,1}', [RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InstructOnsets_run-0',char(string(r)),'.txt'], 'Delimiter', ' ', FileType='text');
           if isempty(Runs{r,1}.onset(Runs{r,1}.corrvsincor==0)');
                    writematrix('*',[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InCorr_run-0',char(string(r)),'.txt'],FileType='text');
                else 
                    writematrix(Runs{r,1}.onset(Runs{r,1}.corrvsincor==0)',[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InCorr_run-0',char(string(r)),'.txt'], 'Delimiter', ' ',FileType='text');
           end
        else dispOut{1,1}(end-8:end) == 'Version 2';
            writematrix(Inst{r,1}', [RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InstructOnsets_run-0',char(string(r)),'.txt'], 'Delimiter', ' ',FileType='text');
            if isempty(Runs{r,1}.onset(Runs{r,1}.corrvsincor==0)');
                    writematrix('*',[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InCorr_run-0',char(string(r)),'.txt'],FileType='text');
                else 
                    writematrix(Runs{r,1}.onset(Runs{r,1}.corrvsincor==0)',[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_InCorr_run-0',char(string(r)),'.txt'], 'Delimiter', ' ',FileType='text');
           end
        end
        clear x y x_str Instruct w i row col isWord j m m2 
        cd (RawDataFolder)
        if r == 2
            allruns = [Runs{1,1};Runs{2,1}];
            writetable(allruns,[RawDataFolder,'working/sub-',subID,'/sub-',subID,'_ses-1_task-siat_allruns_eventss.csv'],"FileType","text", 'Delimiter', ' ');
     end
 
    end
    
  end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %this section is run after you have ran the above section of the script
 %and copies over the tsv files from the source folder to the deriviatves
 %folder
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fMRIprepDone == 1
     if cp_source2deriv == 1;
         for S = 1:size(Subjects,2)
             subID = num2str(Subjects(S)); %subject id
             cd ([SourcedataFolder,subID,'/ses-1/func/'])
             fmriprepoutfold = [DerivativesFolder,fmriprepFolder,subID,'/ses-1/func/'];
             copyfile('*events.tsv', fmriprepoutfold )
             copyfile('*events.csv', fmriprepoutfold)
         end
     end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %this section creates the rest of the afni stimuli files that are needed
 %to conduct a "full model" of the siat.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MakeAFNIstim == 1;
        for S = 1:size(Subjects,2)
            subID = num2str(Subjects(S)); %subject id
            fmriprepoutfold = [DerivativesFolder,fmriprepFolder,subID,'/ses-1/func/'];
            cd(fmriprepoutfold)
            if ~exist([DerivativesFolder,'afni/regressors/sub-',subID,'/'], 'dir');    mkdir([DerivativesFolder,'afni/regressors/sub-',subID,'/']);    end
            cd([DerivativesFolder,'afni/regressors/sub-',subID,'/'])
            try
                movefile([RawDataFolder,'working/sub-',subID,'/','*.csv'],[DerivativesFolder,'afni/regressors/sub-',subID,'/'])
                movefile([RawDataFolder,'working/sub-',subID,'/','*.txt'],[DerivativesFolder,'afni/regressors/sub-',subID,'/'])
            catch
            end
                fileNames = dir (['*events.csv']);
            for r = 1:size(fileNames,1)
                runNum = ['run-0',num2str(r)];
                Run = readtable(fileNames(r).name);
                if RegressorType == "Full_Model"
                        if Run.version(1) == 1
                        Run = readtable(fileNames(r).name);
                        %Practice regressors
                        %DeathWords During Pratice
                        if isempty(Run.onset((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))');
                            if MakeRegressors == 1; writematrix('*',['sub-',subID,'_DeathWords_prac_',runNum]);else DeathWords_prac = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_prac_',runNum], 'Delimiter', ' ');
                            else DeathWords_prac = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1)))];end
                        end
                        %LifeWords During Pratice
                        if isempty(Run.onset((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))');
                            if MakeRegressors == 1; writematrix('*',['sub-',subID,'_LifeWords_prac_',runNum]);else LifeWords_prac = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_prac_',runNum], 'Delimiter', ' ');
                            else LifeWords_prac = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))),std(Run.responsetime((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1)))];end
                        end
                        %MeWords During Pratice
                        if isempty(Run.onset((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_prac_',runNum]);else MeWords_prac = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))',['sub-',subID,'_MeWords_prac_',runNum], 'Delimiter', ' ');
                            else MeWords_prac = [nanmean(Run.responsetime((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&Run.block==2&Run.corrvsincor==1)))];end                        
                        end
                        %OtherWords During Pratice
                        if isempty(Run.onset((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_prac_',runNum]);else OtherWords_prac = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_prac_',runNum], 'Delimiter', ' ');
                            else OtherWords_prac = [nanmean(Run.responsetime((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&Run.block==2&Run.corrvsincor==1)))];end
                        end
                        %Combined blocks
                        %Me and Death blocks or Me and Life blocks
                        %Death words during blocks where Me and Death are together
                        if isempty(Run.onset((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_DeathWords_Me_',runNum]);else DeathWords_Me = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_Me_',runNum], 'Delimiter', ' ');
                            else DeathWords_Me = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))];end
                        end
                        %Life words during blocks where Me and Life are together
                        if isempty(Run.onset((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_LifeWords_Me_',runNum]);else LifeWords_Me = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_Me_',runNum], 'Delimiter', ' ');
                            else LifeWords_Me = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                        end
                                        if isempty(Run.onset((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_Death_',runNum]);else MeWords_Death = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_MeWords_Death_',runNum], 'Delimiter', ' ');
                            else MeWords_Death = [nanmean(Run.responsetime((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))]; end
                        end
                        %Other words during blocks where Me and Life are together
                        if isempty(Run.onset((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_Death_',runNum]);else OtherWords_Death = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_Death_',runNum], 'Delimiter', ' ');
                            else OtherWords_Death = [nanmean(Run.responsetime((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                        end
                        %Other and Death blocks or Other and Life blocks
                        %Death words during blocks where Other and Death are together
                        if isempty(Run.onset((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_DeathWords_Other_',runNum]);else DeathWords_Other = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_Other_',runNum], 'Delimiter', ' ');
                            else DeathWords_Other = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                        end
                        %Life words during blocks where Other and Life are together
                        if isempty(Run.onset((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_LifeWords_Other_',runNum]);else LifeWords_Other = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_Other_',runNum], 'Delimiter', ' ');
                            else LifeWords_Other = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))]; end
                        end
                        %Me words during blocks where Other and Death are together
                        if isempty(Run.onset((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_Life_',runNum]);else MeWords_Life = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_MeWords_Life_',runNum], 'Delimiter', ' ');
                            else MeWords_Life = [nanmean(Run.responsetime((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                        end
                        %Other words during blocks where Other and Life are together
                        if isempty(Run.onset((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                            if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_Life_',runNum]);else OtherWords_Life = [0,0];end
                        else
                            if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_Life_',runNum], 'Delimiter', ' ');
                            else OtherWords_Life = [nanmean(Run.responsetime((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))]; end
                        end
                        Version = 1;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %this section calculates the dscore and reaction time variables for each
 %subject, then creates a combined file for all particpants data. This is
 %for version 1 of the siat
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if CalcDScore ==1
                             
                                A1mean = nanmean(Run.responsetime(Run.block==3));
                                A2mean = nanmean(Run.responsetime(Run.block==4));
                                B1mean = nanmean(Run.responsetime(Run.block==6));
                                B2mean = nanmean(Run.responsetime(Run.block==7));
                                SD4Latency1 = nanstd([Run.responsetime(Run.block==3);Run.responsetime(Run.block==6)]);
                                SD4Latency2 = nanstd([Run.responsetime(Run.block==4);Run.responsetime(Run.block==7)]);
                                Dscore1 = (B1mean-A1mean)/SD4Latency1; 
                                Dscore2 = (B2mean-A2mean)/SD4Latency2;
                                Dscore = mean([Dscore1,Dscore2]);
                                Subj = str2num(subID);
                                OutTable = table(Subj,r,Version,A1mean,A2mean,B1mean,B2mean,SD4Latency1,SD4Latency2,Dscore1,Dscore2,Dscore);
                                writetable(OutTable,['sub-',subID,'_dscore_run',num2str(r),'_Version',num2str(Version),'.csv']);
                                DScores = [DScores;table2array(OutTable)];
                                
                                clear A1mean A2mean B2mean B1mean SD4Latency1 SD4Latency2 Dscore1 Dscore2 Dscore OutTable
                             
                        end
                        if CalcRTS ==1
                            Subj = str2num(subID);
                            RTS=[RTS;str2num(subID),r,Version,DeathWords_prac,LifeWords_prac,MeWords_prac,OtherWords_prac,DeathWords_Me,LifeWords_Me,DeathWords_Other,LifeWords_Other,MeWords_Death,OtherWords_Death,MeWords_Life,OtherWords_Life];
                            OutTableR = table(Subj,r,Version,DeathWords_prac,LifeWords_prac,MeWords_prac,OtherWords_prac,DeathWords_Me,LifeWords_Me,DeathWords_Other,LifeWords_Other,MeWords_Death,OtherWords_Death,MeWords_Life,OtherWords_Life);
                            writetable(OutTableR,['sub-',subID,'_RTS_run',num2str(r),'_Version',num2str(Version),'.csv']);
                        end
             elseif Run.version(1) == 2
                    if isempty(Run.onset((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))');
                        if MakeRegressors == 1;writematrix('*',['sub-',subID,'_DeathWords_prac_',runNum]);else DeathWords_prac = [0,0];end
                    else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_prac_',runNum], 'Delimiter', ' ');
                    else DeathWords_prac = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==1&(Run.block==1|Run.block==5)&Run.corrvsincor==1)))]; end
                     end
                %LifeWords During Pratice
                if isempty(Run.onset((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_LifeWords_prac_',runNum]);else LifeWords_prac = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_prac_',runNum], 'Delimiter', ' ');
                    else LifeWords_prac = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==2&(Run.block==1|Run.block==5)&Run.corrvsincor==1)))]; end
                end
                %MeWords During Pratice
                if isempty(Run.onset((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_prac_',runNum]);else MeWords_prac = [0,0]; end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))',['sub-',subID,'_MeWords_prac_',runNum], 'Delimiter', ' ');
                    else MeWords_prac = [nanmean(Run.responsetime((Run.stimtype==3&Run.block==2&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&Run.block==2&Run.corrvsincor==1)))]; end
                end
                %OtherWords During Pratice
                if isempty(Run.onset((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_prac_',runNum]);else OtherWords_prac = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_prac_',runNum], 'Delimiter', ' ');
                    else OtherWords_prac = [nanmean(Run.responsetime((Run.stimtype==4&Run.block==2&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&Run.block==2&Run.corrvsincor==1)))]; end
                end
                %Combined blocks
                %Me and Death blocks or Me and Life blocks
                %Death words during blocks where Me and Death are together
                if isempty(Run.onset((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_DeathWords_Me_',runNum]);else DeathWords_Me = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_Me_',runNum], 'Delimiter', ' ');
                    else DeathWords_Me = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanmean(Run.responsetime((Run.stimtype==1&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                end
                %Life words during blocks where Me and Life are together
                if isempty(Run.onset((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_LifeWords_Me_',runNum]);else LifeWords_Me = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_Me_',runNum], 'Delimiter', ' ');
                    else LifeWords_Me = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==2&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))]; end
                end
                %Me words during blocks where Me and Death are together
                if isempty(Run.onset((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_Death_',runNum]);else MeWords_Death = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_MeWords_Death_',runNum], 'Delimiter', ' ');
                    else MeWords_Death = [nanmean(Run.responsetime((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                end
                %Other words during blocks where Me and Life are together
                if isempty(Run.onset((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_Death_',runNum]); else OtherWords_Death = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_Death_',runNum], 'Delimiter', ' ');
                    else OtherWords_Death = [nanmean(Run.responsetime((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))];end
                end
                %Other and Death blocks or Other and Life blocks
                %Death words during blocks where Other and Death are together
                if isempty(Run.onset((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                   if MakeRegressors == 1;writematrix('*',['sub-',subID,'_DeathWords_Other_',runNum]);else DeathWords_Other = [0,0];end
                else
                   if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_DeathWords_Other_',runNum], 'Delimiter', ' ');
                   else DeathWords_Other = [nanmean(Run.responsetime((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==1&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))];end
                    
                end
                %Life words during blocks where Other and Life are together
                if isempty(Run.onset((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_LifeWords_Other_',runNum]);else LifeWords_Other = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_LifeWords_Other_',runNum], 'Delimiter', ' ');
                    else LifeWords_Other = [nanmean(Run.responsetime((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==2&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))]; end
                end
                if isempty(Run.onset((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_MeWords_Life_',runNum]);else MeWords_Life = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))',['sub-',subID,'_MeWords_Life_',runNum], 'Delimiter', ' ');
                    else MeWords_Life = [nanmean(Run.responsetime((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==3&(Run.block==3|Run.block==4)&Run.corrvsincor==1)))];end
                    
                end
                %Other words during blocks where Other and Life are together
                if isempty(Run.onset((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))');
                    if MakeRegressors == 1;writematrix('*',['sub-',subID,'_OthersWords_Life_',runNum]);else OtherWords_Life = [0,0];end
                else
                    if MakeRegressors == 1;writematrix(Run.onset((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))',['sub-',subID,'_OthersWords_Life_',runNum], 'Delimiter', ' ');
                    else OtherWords_Life = [nanmean(Run.responsetime((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1))),nanstd(Run.responsetime((Run.stimtype==4&(Run.block==6|Run.block==7)&Run.corrvsincor==1)))];end
                end
                Version = 2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %this section calculates the dscore and reaction time variables for each
 %subject, then creates a combined file for all particpants data. This is
 %for version 2 of the siat
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 if CalcDScore ==1
                             
                                A1mean = nanmean(Run.responsetime(Run.block==6));
                                A2mean = nanmean(Run.responsetime(Run.block==7));
                                B1mean = nanmean(Run.responsetime(Run.block==3));
                                B2mean = nanmean(Run.responsetime(Run.block==4));
                                SD4Latency1 = nanstd([Run.responsetime(Run.block==6);Run.responsetime(Run.block==3)]);
                                SD4Latency2 = nanstd([Run.responsetime(Run.block==7);Run.responsetime(Run.block==4)]);
                                Dscore1 = (B1mean-A1mean)/SD4Latency1; 
                                Dscore2 = (B2mean-A2mean)/SD4Latency2;
                                Dscore = mean([Dscore1,Dscore2]);
                                Subj = str2num(subID);
                                OutTable = table(Subj,r,Version,A1mean,A2mean,B1mean,B2mean,SD4Latency1,SD4Latency2,Dscore1,Dscore2,Dscore);
                                writetable(OutTable,['sub-',subID,'_dscore_run',num2str(r),'_Version',num2str(Version),'.csv']);
                                DScores = [DScores;table2array(OutTable)];
                                
                                clear A1mean A2mean B2mean B1mean SD4Latency1 SD4Latency2 Dscore1 Dscore2 Dscore Subj OutTable
                            
                 end
                 if CalcRTS ==1
                        Subj = str2num(subID);
                            RTS=[RTS;str2num(subID),r,Version,DeathWords_prac,LifeWords_prac,MeWords_prac,OtherWords_prac,DeathWords_Me,LifeWords_Me,DeathWords_Other,LifeWords_Other,MeWords_Death,OtherWords_Death,MeWords_Life,OtherWords_Life];
                            OutTableR = table(Subj,r,Version,DeathWords_prac,LifeWords_prac,MeWords_prac,OtherWords_prac,DeathWords_Me,LifeWords_Me,DeathWords_Other,LifeWords_Other,MeWords_Death,OtherWords_Death,MeWords_Life,OtherWords_Life);
                            writetable(OutTableR,['sub-',subID,'_RTS_run',num2str(r),'_Version',num2str(Version),'.csv']);

                 end
             else
                 NameIssue = 'There is no version lable in file, please check the version and add a varible to the file wtih version information.'
             end
         end
                

                
        end
         if r ==2
             command = ['cat sub-',subID,'_InCorr_run-01.txt sub-',subID,'_InCorr_run-02.txt > sub-',subID,'_InCorr_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_InstructOnsets_run-01.txt sub-',subID,'_InstructOnsets_run-02.txt > sub-',subID,'_InstructOnsets_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_DeathWords_prac_run-01.txt sub-',subID,'_DeathWords_prac_run-02.txt > sub-',subID,'_DeathWords_prac_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_LifeWords_prac_run-01.txt sub-',subID,'_LifeWords_prac_run-02.txt > sub-',subID,'_LifeWords_prac_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_MeWords_prac_run-01.txt sub-',subID,'_MeWords_prac_run-02.txt > sub-',subID,'_MeWords_prac_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_OthersWords_prac_run-01.txt sub-',subID,'_OthersWords_prac_run-02.txt > sub-',subID,'_OthersWords_prac_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_DeathWords_Me_run-01.txt sub-',subID,'_DeathWords_Me_run-02.txt > sub-',subID,'_DeathWords_Me_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_LifeWords_Me_run-01.txt sub-',subID,'_LifeWords_Me_run-02.txt > sub-',subID,'_LifeWords_Me_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_MeWords_Death_run-01.txt sub-',subID,'_MeWords_Death_run-02.txt > sub-',subID,'_MeWords_Death_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_OthersWords_Death_run-01.txt sub-',subID,'_OthersWords_Death_run-02.txt > sub-',subID,'_OthersWords_Death_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_DeathWords_Other_run-01.txt sub-',subID,'_DeathWords_Other_run-02.txt > sub-',subID,'_DeathWords_Other_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_LifeWords_Other_run-01.txt sub-',subID,'_LifeWords_Other_run-02.txt > sub-',subID,'_LifeWords_Other_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_MeWords_Life_run-01.txt sub-',subID,'_MeWords_Life_run-02.txt > sub-',subID,'_MeWords_Life_allRuns.txt'];system(command,'-echo');
             command = ['cat sub-',subID,'_OthersWords_Life_run-01.txt sub-',subID,'_OthersWords_Life_run-02.txt > sub-',subID,'_OthersWords_Life_allRuns.txt'];system(command,'-echo');
         end


    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %this section creates the combined dscore and rts files for all subjects
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if CalcRTS == 1       
 RTS = array2table(RTS);
 RTS.Properties.VariableNames = {'ID','Run','Version','DeathWords_prac_mean','DeathWords_prac_sd',...
     'LifeWords_prac_mean','LifeWords_prac_sd','MeWords_prac_mean','MeWords_prac_sd',...
     'OtherWords_prac_mean','OtherWords_prac_sd','DeathWords_Me_mean','DeathWords_Me_sd',...
     'LifeWords_Me_mean','LifeWords_Me_sd','DeathWords_Other_mean','DeathWords_Other_sd',...
     'LifeWords_Other_mean','LifeWords_Other_sd','MeWords_Death_mean','MeWords_Death_sd',...
     'OtherWords_Death_mean','OtherWords_Death_sd','MeWords_Life_mean','MeWords_Life_sd',...
     'OtherWords_Life_mean','OtherWords_Life_sd'};
  writetable(RTS,[BehavDataFolder,RTS_filename,'_All.csv'])
 writetable(RTS(RTS.Run==1,:),[BehavDataFolder,RTS_filename,'_Run1.csv'])
 else;end;
 if CalcDScore ==1
     DScores = array2table(DScores);
    DScores.Properties.VariableNames = {'Subj','Run','Version','A1mean','A2mean','B1mean','B2mean',...
     'SD4Latency1','SD4Latency2','Dscore1','Dscore2','Dscore'};
     writetable(DScores,[BehavDataFolder,DS_filename,'_All.csv']);
    writetable(DScores(DScores.Run==1,:),[BehavDataFolder,DS_filename,'_Run1.csv'],'WriteVariableNames',true);
 else;end;

     end
 end



  
      
 

