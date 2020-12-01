%%%

%% MARCH BETTER EXP , experimental data I compare with

mainDataDir = '/Volumes/storage/Eleana/ImageData2020/XMAS_12_23_19_esio17_wnts_axin/images';
dDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');
nChannels = 3;
[RealradialAvgNuc] = RealColoniesAnalysis(dDir,nChannels);

%anonymous function 
RealradialAvgNuc;

%trim down tommatch size

NewRealradialAvgNuc = {zeros(6,3),zeros(6,3),zeros(6,3),zeros(6,3),zeros(6,3)};
for k =1:size(RealradialAvgNuc,2)
RealradialAvgNuc{k} = RealradialAvgNuc{k}(1:66, 1:3);

for i = 11:11:66
     NewRealradialAvgNuc{k}(i/11,1:3) = RealradialAvgNuc{k}(i, 1:3);
end

end

% save(['result_data_run_number_' num2str(n) ',mat'],'s')
save(fullfile(dDir,"RealData.mat"), "NewRealradialAvgNuc" )
clear all; clc

mainDataDir = '/Users/elenirea/Documents/WarmflashWNTproject/ImageData/XMAS_12_23_19_esio17_wnts_axin/stiched4images/images';
dDir = mainDataDir;  % fullfile(mainDataDir, 'LSM10x');
load(fullfile(dDir,"RealData.mat"))
NewRealradialAvgNuc;

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


costall = @(x)costFunsolveFun(NewRealradialAvgNuc, x);
lb = [10,0.01*ones(1,7)];
ub = [15,0.3*ones(1,7)];
ga(costall,8, [],[],[],[],lb, ub )

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


8exceptionally well - it worked
 costall([10.0000    0.0100    0.0100    0.0100    0.3000    0.0432    0.3000    0.0100])
