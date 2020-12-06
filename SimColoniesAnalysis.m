function [SimradialAvgNuc]  = SimColoniesAnalysis(dapi) 


fInFormat ='u%d_%dh.tif'; %input
fOutFormat = 'u_%dh.tif'; %output
outDirec = '/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling/outPutODE45_2020/images';
inDirec = '/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling/outPutODE45_2020/';
% IWP2 40:4:252
% XMAS 4:4:996


img = zeros(100,100,3);  
% img(:,:,1) = dapi ;

for img_num = [2,4,6,8]

    for jj = 1:2
        fname = fullfile(inDirec,sprintf(fInFormat,jj,img_num));
        if exist(fname,'file')
            img(:,:,jj) = imread(fname);
        end
    end
        img(:,:,3) = dapi ;
    a = mean(mean(img(:,:,2)));
    img(:,:,3)  = img(:,:,3).*a;
    
    if exist('img','var') && size(img,3) == 3
        fOut = fullfile(outDirec,sprintf(fOutFormat,img_num));
        imwrite(img,fOut);
    end
    
    clear img
end
%%

warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

mainDataDir = '/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling/outPutODE45_2020/images';
dataDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');

fInFormat = 'u_%dh.tif';
filename = sprintf(fInFormat, 1);
colRadius = 30;
fileName = [dataDir filesep filename];
meta = MetadataMicropattern(fileName);
meta.colMargin = 1;
% fix some stuff manually
meta.xres = 0.9;
meta.yres = 0.9;
meta.nChannels = 3;
meta.nZslices = 1;
meta.nTime = 1;
% meta.excitationWavelength = wavs

meta.colRadiiMicron = colRadius;
meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);

fileRange = folderFilesFromKeyword(dataDir,'u_');

% u2
filenrs = {fileRange(4), fileRange(3),fileRange(2), fileRange(1) }; %get actual file numbers

meta.channelLabel = { 'WNT', 'AXIN2','DAPI'};
% filenrs = {4:5 9:11 15:17 21:23 27:29 33:35};
% meta.channelLabel = {'DAPI','SIX1','PAX3','SNAI1'};

meta.conditions = ["48h", "30h", "24h", "12h"];

save(fullfile(dataDir,'metaData.mat'),'meta');

DAPIChannel = 3;
clear colonies
%% process (save DAPI channel MIP for segmentation)

% although findColonies can find multiple colonies in a single image,
% for the LSM there is one per image, which this code below assumes
findColParam = struct('sclose', 6, 'sopen', 8, 'checkcontained', false,...
    'minArea', [],'convhull', true);
% close all;
colonies(numel([filenrs{:}])) = Colony;

for coli = [filenrs{:}]
    
    % cleanScale in micron
    param = {'DAPIChannel',DAPIChannel, 'colID',coli, 'adjustmentFactor',0.5,'clparameters',findColParam};
    filename = sprintf(fInFormat, coli); %i want to see the name
    colony = processOneSimColonyImage(filename, dataDir, param);
    colonies(coli) = colony;
    colonies(coli).setID(coli);
    % setID overwrites filename, which was the right thing for epi, not
    % here
    colonies(coli).filename = filename;
end

save(fullfile(dataDir,'colonies'), 'colonies');

%% show plot of different conditions side by side

load(fullfile(dataDir,'colonies'), 'colonies');

% first normalize by DAPI, then scale all profiles from 0 to 1

for i = 1:numel(meta.conditions)
    coloniesCombined{i} = colonies(filenrs{i});
end

%inside function to make averages and plot

doubleNormalize = 0;
n = numel(coloniesCombined);
m = 1;
radialAvgNuc = {};
r = {};
minI = Inf*(1:meta.nChannels);
maxI = 0*(1:meta.nChannels);

for i = 1:n
    radialAvg = makeAveragesNoSegmentation(...
                    meta, colRadius, DAPIChannel, coloniesCombined{i});
%     chans = 1:length(meta.channelLabel);
%     chansToPlot = setdiff(chans,DAPIChannel);
    chansToPlot = 1;
    
    if ~isempty(DAPIChannel)
        radialAvgNuc{i} = radialAvg.nucAvgDAPINormalized;
    else
        radialAvgNuc{i} = radialAvg.nucAvg;
    end
    r{i} = radialAvg.r;    
    
    % for overall normalization
    % throw out 2 bins from edge when setting LUT
    % to prevent setting minimum by areas without cells
    Imargin = 0; 
    minI = min(minI, min(radialAvgNuc{i}(1:end-Imargin,:)));
    maxI = max(maxI, max(radialAvgNuc{i}(1:end-Imargin,:)));
end

if doubleNormalize == 1
    for i = 1:n
        for ci = 1:meta.nChannels
            radialAvgNuc{i}(:,ci) = (radialAvgNuc{i}(:,ci) - minI(ci))/(maxI(ci)-minI(ci)); %for DAPI WE DO NOT CARE - 1-1 =0
        end
    end
end


figure('Position',[0 0 1600 300]);
for i = 1:n
    subplot_tight(m,n,i,0.02)
    plot(r{i}, radialAvgNuc{i}(:,chansToPlot),'.-','LineWidth',3)
    axis([min(r{i}) max(r{i}) 0 1]);
    legend(meta.channelLabel(chansToPlot));
    title(meta.conditions{i})
    axis square
    if i > 1
        legend off;
    end
end
    saveas(gcf,fullfile(dataDir,'radialProfilesSimCombinedAXIN2.png'));
%     radialAvgNuc{1}
%     radialAvgNuc{4}
    SimradialAvgNuc = radialAvgNuc;
end