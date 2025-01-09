clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Paramaters to run anlaysis, modify to fit needs                  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
studyDir = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI';
atlasDir = '/quobyte/nerve-data/dsalat/2389_trt_main/arickels/Experiments/RP010_siat_tfMRI/derivatives/afni/templates/';
contrastDir = '/derivatives/afni/SecondLevelModels/FullMod_Blur6_Polort5_3to44_5SubsRem/LME_wordBYcong/';
clusterMaskName = 'Clust_mask_CongbyWordF_p1neg3_clusterize38+tlrc';
FLDir = '/derivatives/afni/FirstLevelModels/FullMod_Blur6_Polort5/';
atlasName = '7net200WithAmygHip';
output_name = 'WholeBrain_coef_CongBYwordF_p1neg3_clust38'
ROI_indx = [168];%168; 201
ROI_lable = ['RightAmyg'];
output_data = 1; %do you want to output the data into a file
calcYeo = 1; %1 means calculate percentage in each yeo network

%set up what will be extracted and from where
CoEF_TTest = 0; %0 = use coeff, 1 = use tstat
ClustOrWholeMap =  0; %cluster == 1 and whole map == 0;
AtlasOrClustOut = 0; %either calcuate from rois in an atals or from cluster
                     %output from afni
Cont_TorV = 1; %are you going to use the coef and t stat from either the ttest 
               % done in the script (t = 0), or you plan to look at each
               % regressor (v = 1).
doBar = 0; %this indicates if you want to do bar graphs with output 

SingleCont = []; %identify the indxe of the exact index you want to do from ContIndx
%set up subject info
removeSubs = [3,25,31,36,40];
Subjects = [3:44];

%clean up cluster outupt
removeCluster = 0; %1 == removeing clusters, 0 = keep all clusteres
RM_clustNum = [0]; %cluster numbers to remove

%list of all contrast names and indexes for FullModel
ClustDir{1,1} = [studyDir,contrastDir];
SubDir{1,1} = [studyDir,FLDir];
%the ContName variable is useful when I was doing indivudal t-tests
if Cont_TorV == 0
    if CoEF_TTest == 0
        ContNames{1,1} = [{'DeathWordsIncogr-DeathWordsCongr_coef'},...
                {'LifeWordsIncogr-LifeWordsCongr_coef'},...
                {'MeWordsIncogr-MeWordsCongr_coef'},...
                {'NotMeWordsIncogr-NotMeWordsCongr_coef'},...
                {'Incong_VS_Cong_coef'}];
        ContIndx = [44,47,50,53,56]; %coeffecents
    else
        ContNames{1,1} = [{'DeathWordsIncogr-DeathWordsCongr_tstat'},...
                {'LifeWordsIncogr-LifeWordsCongr_tstat'},...
                {'MeWordsIncogr-MeWordsCongr_tstat'},...
                {'NotMeWordsIncogr-NotMeWordsCongr_tstat'},...
                {'Incong_VS_Cong_tstat'}];
        ContIndx = [45,48,51,54,57]; %tstat
    end
else 
    if CoEF_TTest ==0
        ContNames{1,1} = [{'Incogr-Cong_coef'},{'DeathWordsIncogr_coef'},{'DeathWordsCongr_coef'},...
                {'LifeWordsIncogr_coef'},{'LifeWordsCongr_coef'},{'MeWordsIncogr_coef'},{'MeWordsCongr_coef'},...
                {'OtherWordsIncogr_coef'},{'OtherWordsCongr_coef'}];
        ContIndx = [56,5,8,17,14,26,23,32,35]; %coeffecents
    else
        ContNames{1,1} = [{'Incogr-Cong_tstat'},{'DeathWordsIncogr_tstat'},{'DeathWordsCongr_tstat'},...
                {'LifeWordsIncogr_tstat'},{'LifeWordsCongr_tstat'},{'MeWordsIncogr_tstat'},{'MeWordsCongr_tstat'},...
                {'OtherWordsIncogr_tstat'},{'OtherWordsCongr_tstat'}];
        ContIndx = [57,6,9,18,15,27,24,33,36]; %tstat
    end
    %%% REMEMBER WHEN IDENTIFYING SUBBRICKS THAT THE NUMBERING STARTS AT 0
    %%% BUT MATLAB STARTS AT 1. SO IF THE SUBBRICK IS 46, WHEN YOU
    %%% IDENTIFY IT IN MATLAB IT WOULD BE 47.
end
WordLength = [3,9,6,4;7,7,4,5;8,4,2,4;7,6,1,5;8,5,4,4];
MeanWordLength = mean(WordLength);
if ~isempty(SingleCont)
    ContIndx = ContIndx(SingleCont);
else
end
%Load Data File and Calculate additional variables.
DC = readtable([studyDir,'/backup_data/FinalDataFile.csv']);
    DC.DIncMinDCong = DC.DeathWords_Me_mean-DC.DeathWords_Other_mean;
    DC.LIncMinLCong = DC.LifeWords_Other_mean-DC.LifeWords_Me_mean;
    DC.MIncMinMCong = DC.MeWords_Death_mean-DC.MeWords_Life_mean;
    DC.OIncMinOCong = DC.OtherWords_Life_mean-DC.OtherWords_Death_mean;
    DC.DS_grouped = DC.Dscore>0;

RS = removeSubs-2;
Subjects(RS)=[];
DC(RS,:) = [];
%%%%Load cluster mask
[~,mask,info] = BrikLoad ([ClustDir{1,1},clusterMaskName]);
for i = 1:size(mask,3)
    mask(:,:,i)=flipud(mask(:,:,i));
end
NumberOfClusters = max(unique(mask));
NumberOfClusters
clusterMaskName

if AtlasOrClustOut == 1;
    %atlas =load([atlasDir,atlasName]);
    %atlas2 = imresize3(atlas.UseMask,[97,115,97],"Method","nearest");
    %mask = atlas2;
    holdit = load([atlasDir,atlasName]);
    mask = holdit.atlas; clear holdit
   NumberOfClusters = max(unique(mask));
    NumberOfClusters
    
end

%%%remove clusters we dont want
if removeCluster ==1
    for i = 1:length(RM_clustNum)
        mask(mask==RM_clustNum(i)) = 0;
    end
end
clustNum = unique(mask);clustNum(1) = [];

%calc percentage in yeo 7 network
if calcYeo == 1
    holdit = load([atlasDir,atlasName]);
    Atlas = holdit.atlas; clear holdit
    atlas = reshape(Atlas,size(Atlas,1)*size(Atlas,2),size(Atlas,3));
    for i = [1:14,101:115];    atlas(atlas==i) = 1;    end
    for i = [15:30,116:134];    atlas(atlas==i) = 2;    end
    for i = [31:43,135:147];    atlas(atlas==i) = 3;    end
    for i = [44:54,148:158,201:204];    atlas(atlas==i) = 4;    end
    for i = [55:60,159:164];    atlas(atlas==i) = 5;    end
    for i = [61:73,165:181];    atlas(atlas==i) = 6;    end
    for i = [74:100,182:200];    atlas(atlas==i) = 7;    end
    Atlas = reshape(atlas,size(Atlas,1),size(Atlas,2),size(Atlas,3));


    if ClustOrWholeMap == 0 
        mask2 = mask>0;
        test = mask2.*Atlas;
        Visual = ((sum(sum(sum(test==1))))/(sum(sum(sum(mask>0)))))*100
        SM = ((sum(sum(sum(test==2))))/(sum(sum(sum(mask>0)))))*100
        DAN = ((sum(sum(sum(test==3))))/(sum(sum(sum(mask>0)))))*100
        VAN = ((sum(sum(sum(test==4))))/(sum(sum(sum(mask>0)))))*100
        Limb = ((sum(sum(sum(test==5))))/(sum(sum(sum(mask>0)))))*100
        CCN = ((sum(sum(sum(test==6))))/(sum(sum(sum(mask>0)))))*100
        DMN = ((sum(sum(sum(test==7))))/(sum(sum(sum(mask>0)))))*100
    else
        for x = 1:(max(unique(mask)))
            mask2 = mask==x;
            test = mask2.*Atlas;
            clustA(x,1) = ((sum(sum(sum(test==1))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,2) = ((sum(sum(sum(test==2))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,3) = ((sum(sum(sum(test==3))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,4) = ((sum(sum(sum(test==4))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,5) = ((sum(sum(sum(test==5))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,6) = ((sum(sum(sum(test==6))))/(sum(sum(sum(mask2>0)))))*100;
            clustA(x,7) = ((sum(sum(sum(test==7))))/(sum(sum(sum(mask2>0)))))*100;
        end
    end
else
end


%%Apply mask to individual subject stats file
if size(clustNum) >0
    for j = 1:length(Subjects)
        SubIndx = strcat("Particpant ", num2str(j))
        [~,subBrain,info] = BrikLoad([SubDir{1,1},'sub-',num2str(Subjects(j)),'/stats.sub-',num2str(Subjects(j)),'_run-01+tlrc']);
        subBrain(subBrain==0)=nan;
        subBrain(isinf(subBrain))=nan;
        
         for i = 1:length(ContIndx)
             %%if you are only running one contrast, this sets it to the
             %%correct index
                if ~isempty(SingleCont)
                    i = SingleCont
                else
                end
                
                ContNames{1,1}{i}
                if AtlasOrClustOut==1
                    subBrain_contrast_mean = subBrain(:,:,:,ContIndx(i));
                    maskBrainVals=subBrain_contrast_mean(mask==ROI_indx);
                    MaskMean((j),i)=nanmean(maskBrainVals(:));
                    MaskSD((j),i)=nanstd(maskBrainVals(:));
                else
                
                if ClustOrWholeMap == 0; %if you want to do a whole brain average
                            subBrain_contrast_mean = subBrain(:,:,:,ContIndx(i));
                            maskBrainVals=subBrain_contrast_mean(mask>0);
                            MaskMean((j),i)=nanmean(maskBrainVals(:));
                            MaskSD((j),i)=nanstd(maskBrainVals(:));
                            


        
                else %if you want to do this anlaysis over every cluster
                            subBrain_contrast_mean = subBrain(:,:,:,ContIndx(i));
                            for x = 1:(max(unique(mask)))
                                maskBrainVals=subBrain_contrast_mean(mask==x);
                                MaskMean((j),i,x)=nanmean(maskBrainVals(:));
                                MaskSD((j),i,x)=nanstd(maskBrainVals(:));
                            end
                       
                end
                end
         end
    end
else   
end

if doBar == 1
%%additional stats
try
%bargraphs
barMeans = mean(MaskMean(:,:,1)); 
bar(barMeans(2:end))
for i = 1:13
subplot(3,5,i)
barMeans = mean(MaskMean(:,:,i)); 
bar(barMeans(2:end))
end

for i = 1:3
subplot(1,3,i)
barMeans = [mean(mean(MaskMean(:,2:3,i),2)),mean(mean(MaskMean(:,4:5,i),2)),...
    mean(mean(MaskMean(:,6:7,i),2)),mean(mean(MaskMean(:,8:9,i),2))];
bar(barMeans)
end

for i = 1:4
subplot(2,2,i)
barMeans = [mean(MaskMean(:,[2,4,6,8],i),2),mean(MaskMean(:,[3,5,7,9],i),2)];
dscore = DC.Dscore;
scatter(barMeans(:,1),dscore, "red","*"); hold on;
scatter(barMeans(:,2),dscore, "blue"); 
h = lsline; 
set(h(1,2),"color","red")
set(h(1,1),"color","blue")
hold off;
[r,p] = corr(barMeans(:,1),dscore,"type","Pearson");
cor_rps(i,1) = r;
cor_rps(i,2) = p;
[r,p] = corr(barMeans(:,2),dscore,"type","Pearson");
cor_rps(i,3) = r;
cor_rps(i,4) = p;
end
catch
end


RTSbyWord  = [mean([DC.DeathWords_Me_mean,DC.DeathWords_Other_mean],2),...
    mean([DC.LifeWords_Me_mean,DC.LifeWords_Other_mean],2),...
    mean([DC.MeWords_Death_mean,DC.MeWords_Life_mean],2),...
    mean([DC.OtherWords_Death_mean,DC.OtherWords_Life_mean],2)];

RTSbyWord_Cong = mean([DC.DeathWords_Other_mean,DC.LifeWords_Me_mean,DC.MeWords_Life_mean,DC.OtherWords_Death_mean]);
RTSbyWord_Incong = mean([DC.DeathWords_Me_mean,DC.LifeWords_Other_mean,DC.MeWords_Death_mean,DC.OtherWords_Life_mean]);
xx = [RTSbyWord_Cong,RTSbyWord_Incong];
ss = [(std(DC.DeathWords_Other_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.LifeWords_Me_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.MeWords_Life_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.OtherWords_Death_mean)/sqrt(length(DC.DeathWords_Other_mean)));...
    (std(DC.DeathWords_Me_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.LifeWords_Other_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.MeWords_Death_mean)/sqrt(length(DC.DeathWords_Other_mean))),...
    (std(DC.OtherWords_Life_mean)/sqrt(length(DC.DeathWords_Other_mean)))]';
 errorbar(xx',ss'); xlim([.80,2.2]) 
else
end
if ClustOrWholeMap == 0
    varNames = [{'CongEff'},{'DW_I'},{'DW_C'},{'LW_I'},{'LW_C'},...
    {'MW_I'},{'MW_C'},{'NW_I'},{'NW_C'},];
    MaskMeanT = array2table(MaskMean,VariableNames=varNames);
else
    NumOfVars = size(ContNames{1,1},2);
    varNames = [{'CongEff_1'},{'DW_I_1'},{'DW_C_1'},{'LW_I_1'},{'LW_C_1'},...
    {'MW_I_1'},{'MW_C_1'},{'NW_I_1'},{'NW_C_1'}];
    mm_temp = MaskMean(:,1:NumOfVars);
    lowDex = 1; upDex = NumOfVars;
    for i = 2:(size(clustNum,1))
        varNames = [varNames,{['CongEff_',num2str(i)]},{['DW_I_',num2str(i)]},...
            {['DW_C_',num2str(i)]},{['LW_I_',num2str(i)]},{['LW_C_',num2str(i)]},...
        {['MW_I_',num2str(i)]},{['MW_C_',num2str(i)]},{['NW_I_',num2str(i)]},...
        {['NW_C_',num2str(i)]}];
        lowDex = lowDex+NumOfVars; upDex = upDex+NumOfVars;
        mm_temp = [mm_temp,MaskMean(:, lowDex:upDex)];
    end
    MaskMeanT = array2table(mm_temp,VariableNames=varNames);
end

if output_data == 1
    try
        writematrix(MaskMeanT,['Output/',output_name,'mean.csv'])
    catch

        writetable(MaskMeanT,['Output/',output_name,'mean_tbl.csv'])
    end
end