
function [RealradialAvgNuc] = RealColoniesAnalysis(dataDir,nChannels)

warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

fInFormat = 'Image%.4d_01.tif';
filename = sprintf(fInFormat, 1);
colRadius = 350;
fileName = [dataDir filesep filename];
meta = MetadataMicropattern(fileName);

% fix some stuff manually
meta.xres = 0.625;
meta.yres = 0.625;
meta.nChannels = nChannels;
meta.nZslices = 1;
meta.nTime = 1;
% meta.excitationWavelength = wavs

meta.colRadiiMicron = colRadius;
meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);

fileRange = folderFilesFromKeyword(dataDir,'Image');

% BLOCK5B = {{1 5,22 26}, {6 10,27 32}, {186 194,33 38} ,{11 16,39 44}, {17 21, 45 50}};  %Block array
% filenrs = {fileRange([1,3:4 ,22:23,26]), fileRange([28:32]), fileRange([186:192,194,33:38]), fileRange([14,16,39:42,44]), fileRange([17:19, 45,46:48:50])}; %get actual file numbers

%BLOCK6 = {{1 5} {11 14} {19 22} {28 32}};  %Block array / removed a bad photo N0. 15
%BLOCK6 = {{51 56,82 87}, {57 64,88 97}, {65 70, 98 104} ,{71 75,105 109}, {76 81, 110 116}};  %Block array
% filenrs = {fileRange([51:56, 82:87]), fileRange([57:64,88:90,92:95,97]), fileRange([69, 99:100, 102:104]), fileRange([71:75,106,109]), fileRange([76:77,79:81, 111:112,114,116]) }; %get actual file numbers

% BLOCK8A = {{117 122,150 155}, {123 130,156 168}, {131 138,169 178} ,{139 144,179 185}, {145 149,193 198}};  %Block array
% filenrs = {fileRange([117,119:121,151:155]), fileRange([124:130,156:162,164:168]), fileRange([131:138,169:178]), fileRange([139:142,144,179:185]), fileRange([149,194,196:198])}; %get actual file numbers

% BLOCKW3 = {{22 26, 82 87, 150 155}, {27 32, 88 97, 156 168}, {33 38, 98 104,169 178} ,{39 44, 105 109, 179 185}, {45 50, 110 116,193 198}};  %Block array;  %Block array
% filenrs = {fileRange([22,23,26, 82:87, 151:155]), fileRange([28:32, 88:90,92:95,97, 156:162,164:168]), fileRange([33:38, 100,102:104,169:178]), fileRange([39:42,44, 106,109, 179:185]), fileRange([45,46,48,49, 111,112, 114,116,194,196:198])}; %get actual file numbers

% BLOCKA2 = {{1 5, 51 56, 117 123}, {6 10,57 64, 123 130}, {186 194, 65 70, 131 138} ,{11 16,71 75, 139 144}, {17 21, 76 81, 145 149}};  %Block array
% fileRange([1,3:4, 51:56, 117,119:121]),
filenrs = { fileRange([57:64, 124:130]), fileRange([186:192,194, 69, 131:136]), fileRange([14,16,72:75, 139:142,144]), fileRange([17:19, 76,77,79:81,149])}; %get actual file numbers

meta.channelLabel = {'DAPI','WNT','AXIN2'};
% filenrs = {4:5 9:11 15:17 21:23 27:29 33:35};
% meta.channelLabel = {'DAPI','SIX1','PAX3','SNAI1'};

meta.conditions = ["48h", "30h", "24h", "12h"];

save(fullfile(dataDir,'metaData.mat'),'meta');

DAPIChannel = 1;

%% process (save DAPI channel MIP for segmentation)

% although findColonies can find multiple colonies in a single image,
% for the LSM there is one per image, which this code below assumes
findColParam = struct('sclose', 6, 'sopen', 8, 'checkcontained', false,...
    'minArea', [],'convhull', true);
% close all;
colonies(numel([filenrs{:}])) = Colony;


for coli = [filenrs{:}]
    
    % cleanScale in micron
    param = {'DAPIChannel',DAPIChannel, 'colID',coli, 'adjustmentFactor', 0.01,'clparameters',findColParam};
    filename = sprintf(fInFormat, coli);
    colony = processOneColonyImage(filename, dataDir, param);
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
doubleNormalize = true;

for i = 1:numel(meta.conditions)
    coloniesCombined{i} = colonies(filenrs{i});
end

coloniesCombined;
%change that to return the averages of all the colonies I input - now only
%coloniesCombined{1}
figure('Position',[0 0 1600 300]);
[radialAvgNuc, r] = plotMultipleAveragesNoSegmentation(meta, colRadius, DAPIChannel,...
    coloniesCombined, meta.conditions, doubleNormalize)
saveas(gcf,fullfile(dataDir,'radialProfilesCombinedAXIN2.png'));

RealradialAvgNuc = radialAvgNuc
end