%%%

%% MARCH BETTER EXP , experimental data I compare with

mainDataDir = '/Volumes/storage/Eleana/ImageData2020/XMAS_12_23_19_esio17_wnts_axin/images';

% mainDataDir = '/Volumes/storage/Eleana/ImageData2020/XMAS_12_23_19_esio17_wnts_axin/XMASclean/IMAGE';
dataDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');

fInFormat = 'Image%.4d_01.tif';
filename = sprintf(fInFormat, 1);
colRadius = 350; % 350 um
fileName = [dataDir filesep filename];
meta = MetadataMicropattern(fileName);
meta.conditions = ["48HC", "48H", "30H","24H", "12H"];
meta.xres =  0.66; %um/p
meta.yres =  0.66;
meta.nChannels = 3;
meta.nZslices = 1;
meta.nTime = 1;
% meta.excitationWavelength = wavs

meta.colRadiiMicron = colRadius;
meta.colRadiiPixel = round(meta.colRadiiMicron/meta.xres);

%%extras
lims = {[0.003 0.02],[0.005 0.04]};
thres = {0.2,0.15, 0.12};
MAINFILENAME = 'WNTsComboXMAS';
chans = {{'WNT5B','WNT3';'WNT5B','WNT3';'WNT5B','WNT3';'WNT5B','WNT3';'WNT5B','WNT3'},{'WNT6','WNT3';'WNT6','WNT3';'WNT6','WNT3';'WNT6','WNT3';'WNT6','WNT3'},{'WNT8A','WNT3';'WNT8A','WNT3';'WNT8A','WNT3';'WNT8A','WNT3';'WNT8A','WNT3'}};
DAPIChannel = 3;


fileRange = folderFilesFromKeyword(dataDir,'Image');

%XMAS
% % BLOCK5B = {{1 5,22 26}, {6 10,27 32}, {186 194,33 38} ,{11 16,39 44}, {17 21, 45 50}};  %Block array
filenrs_WNT5B = {fileRange([2,3:4 ,22:23,26]), fileRange(28:32), fileRange([186:192,33:38]), fileRange([14,16,39:42,44]), fileRange([17:19, 45,46:48:50])}; %get actual file numbers
%
% %BLOCK6 = {{1 5} {11 14} {19 22} {28 32}};  %Block array / removed a bad photo N0. 15
% %BLOCK6 = {{51 56,82 87}, {57 64,88 97}, {65 70, 98 104} ,{71 75,105 109}, {76 81, 110 116}};  %Block array
filenrs_WNT6 = {fileRange([51:56, 82:87]), fileRange([57,58,60:64,88:90,92:95,97]), fileRange([69, 99:100, 102:104]), fileRange([71:75,106,109]), fileRange([76:77,79:81, 111:112,114,116]) }; %get actual file numbers
%
% % BLOCK8A = {{117 122,150 155}, {123 130,156 168}, {131 138,169 178} ,{139 144,179 185}, {145 149,193 198}};  %Block array
filenrs_WNT8A = {fileRange([117,119:121,151:155]), fileRange([124:130,156:161,164:168]), fileRange([131,132,134:136,169:174,176:178]), fileRange([139:142,144,179:185]), fileRange([149])}; %,194,196:198get actual file numbers

% FILLs = {filenrs_WNT5B};FILLsnames = {'filenrs_WNT5B'};
% FILLs = {filenrs_WNT6};FILLsnames = {'filenrs_WNT6'};
% FILLs = {filenrs_WNT8A};FILLsnames = {'filenrs_WNT8A'};
FILLs = {filenrs_WNT5B,filenrs_WNT6, filenrs_WNT8A};FILLsnames = {'filenrs_WNT5B','filenrs_WNT6', 'filenrs_WNT8A'};

[RealradialAvgNuc] = ModRealColoniesAnalysis(meta,dataDir,fInFormat, MAINFILENAME, thres,FILLs,chans,DAPIChannel)

%anonymous function 
RealradialAvgNuc;

%trim down tommatch size
% outdir = '/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling'
% RealradialAvgNuc = load(fullfile(outdir,"RealData.mat"))
a = {zeros(8,3),zeros(8,3),zeros(8,3),zeros(8,3),zeros(8,3)}
NewRealradialAvgNuc = {a a a};

for k = 1:size(NewRealradialAvgNuc,2)
    for v = 1:size(NewRealradialAvgNuc{1},2)
        RRealradialAvgNuc{k}{v} =  RealradialAvgNuc{k}{v}.nucAvg(1:17, 1:3);
        for i = 2:2:17
            NewRealradialAvgNuc{k}{v}(i/2,1:3) = RRealradialAvgNuc{k}{v}(i, 1:3);
        end
    end
end

% save(['result_data_run_number_' num2str(n) ',mat'],'s')

outdir = '/Volumes/storage/Eleana/modelling_gastruloids/NODALmodellling'
save(fullfile(outdir,"RealData.mat"), "NewRealradialAvgNuc" )
clear all; clc
outdir = '/Volumes/storage/Eleana/modelling_gastruloids/NODALmodellling'
load(fullfile(outdir,"RealData.mat"))
% mainDataDir = '/Users/elenirea/Documents/WarmflashWNTproject/ImageData/XMAS_12_23_19_esio17_wnts_axin/stiched4images/images';
% dDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');

%% Run genetic algorithm
% costall([10,0.1,0.1,0.1,0.3,0.001,0.001,0.1])

%9 parameters
% Da = param.Da; 
% ba = param.ba; %basal production A 0.02
% bi = param.bi; %basal production I 0.03
% ka = param.ka; %autoproduction A 0.04
% kb = param.kb; % production due to BMP4 0.2
% kia = param.kia; %inhibition due to I 0.05
% kda = param.kda; % basal degradation A 0.015
% kaa = param.kaa; %production of I 0.03
% kdi = param.kdi; % basal degradation I 0.015
% B = 3; % concentration BMP4

saveInPath = '/Volumes/storage/Eleana/modelling_gastruloids/NODALmodellling/outPutODE45_2020/';

costall = @(x)costFunSolver3D(saveInPath,NewRealradialAvgNuc{1}, x);
lb = [0.01,0.5,0.01*ones(1,5),2,0.01*ones(1,3),5,1]; % there are 14 parameters 2K11_2 = 1 3n11 = 5 8n12 = 5 12K13_2 = 12 13n13 = 5
ub = [0.3,2,0.3*ones(1,5),8,0.01*ones(1,3),15,8];
ga(costall,14, [],[],[],[],lb, ub);
%%%%

%% Steady state values
% f1 = @(A,I) (ba +(A))./((1+kia*I).*(1+ka*A)) - kda*A +Da*del2(A);
% f2 = @(A,I) (kaa*A - kdi*I + bi);
% 
% fff =@(a) [f1(a(1),a(2)) ;f2(a(1),a(2))];
% [equil] = fsolve(fff,[ic(1) ,ic(2)]);
% equil;
%%

%1. Optimization terminated: average change in the fitness value less than options.FunctionTolerance.
% lb = [5,0.001*ones(1,7)];
% ub = [20,0.5*ones(1,7)];
% 19.4509    0.1000    0.3577    0.0230    0.1168    0.4473    0.2200 0.3356


%2.  "cost is"    "0.17277"
% costall([12,0.2,0.1,0.01,0.02,0.1,0.1,1.1])

%3.Optimization terminated: average change in the fitness value less than options.FunctionTolerance. 
% lb = [10,0.001*ones(1,6),0.1];
% ub = [15,0.5*ones(1,6), 1.2];

%4."cost is"    "0.04602"
% 14.9244, 0.0036,  0.4035, 0.1339,  0.1652,  0.4095, 0.4038, 0.8346

% %5."cost is"    "0.090195"
% Optimization terminated: average change in the fitness value less than options.FunctionTolerance.
% 14.9244  0.0038  0.4035  0.1337  0.1652  0.4095  0.4037  0.8346

%6 0.1670 - it is actually good! costall([14,0.2,0.1,0.01,0.16,0.4,0.4,1.1])

%
% Now change both boundaries lb = [10,0.001*ones(1,6),0.1];ub = [15,0.3*ones(1,6), 1.2];
% and cost function

%7 "cost is"    "0.075317"
% costall([14.1744,0.1392,0.1494,0.0504,0.2027,0.1553,0.2242 ,0.6098])
% 14.1744    0.1392    0.1494    0.0504    0.2027    0.1553  0.2242    0.6098


%exceptionally well - it worked
 costall([10.0000    0.0100    0.0100    0.0100    0.3000    0.0432    0.3000    0.0100])
