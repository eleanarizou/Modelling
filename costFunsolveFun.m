
function [costF] = costFunsolveFun(saveInPath, NewRealradialAvgNuc,vec)
disp(vec)
%NewRealradialAvgNuc,
%% Modeling part
radius_colony = 600; %%Radius of colony
L = radius_colony*1.7067;  %domain size
N = 100;
x =linspace(-L,L,N);% (L/N)*(1:N)'; %find the coordinates for an image size of 512
y = x;
[X,Y] = ndgrid(x,x);

R = sqrt(X.^2 + Y.^2);
chi = ((1-tanh((R - radius_colony)/0.1))/2); %whatever positive plot , trash the rest , thresholding

%activators
edgewidth = 100;
edge = (1+tanh((R - radius_colony + edgewidth)/0.1))/2.*(1-tanh((R - radius_colony )/0.1));%.*(1+tanh((R - rad_col)/0.1));
media = ~chi;
%inhibitor
% edgewidth = 200;
% ed = (1+tanh((R - radius_colony + edgewidth)/0.1))/2.*(1-tanh((R - radius_colony )/0.1));%.*(1+tanh((R - rad_col)/0.1));
% 
ic = [0.5, 0.2];
matA = ic(1).*edge;
matI = ic(2).*edge;

%9 parameters
% Da = param.Da;
% ba = param.ba; %basal production A
% bi = param.bi; %basal production I
% ka = param.ka; %autoproduction A
% kia = param.kia; %inhibition due to I
% kda = param.kda; % basal degradation A
% kdi = param.kdi; % basal degradation I
% kaa = param.kaa; %production of I
% kb = param.kb; % production due to BMP4

param.Da = vec(1);
param.ba = vec(2);
param.bi = vec(3);
param.ka = vec(4);
param.kia = vec(5);
param.kda = vec(6);
param.kdi = vec(7);
param.kaa = vec(8);
% param.kb = vec(9);

for i = 2:2:8
    
    f = @(t,y) parametersfun(y, param);
    
    tspan = [i-2;i];
    sol = ode45(f,tspan,[matA(:);matI(:)]);
    y = sol.y(:,end);
    n = length(y)/2;
    A = y(1:n);
    I = y((n+1):end);
    matA = reshape(A,[],sqrt(n));
    matI = reshape(I,[],sqrt(n));
    t = round(sol.x(end));
    imshow(matI,[])
    pause(1)
    imwrite(matA./65,[saveInPath,'u1_',num2str(t),'h.tif']);
    imwrite(matI./65,[saveInPath,'u2_',num2str(t),'h.tif']);
    
end

%% FROM THE SIMULATIONS
dapi = imread('/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling/outPutODE45_2020/dapi0.tif');
dapi = dapi > 30;
dapi = imfill(dapi,'holes');

[SimradialAvgNuc]  = SimColoniesAnalysis(dapi);
SimradialAvgNuc

%% Compute difference and the Cost
% dif = {zeros(6,3),zeros(6,3), zeros(6,3), zeros(6,3)};
% diff = [];
% 
% for y = 1:4
%     dif{y} = abs((NewRealradialAvgNuc{y} - SimradialAvgNuc{y}));
%     diff(:,y) = dif{y}(:,3).^2;
% end
% costF = mean(mean(diff));

%2nd version - compare only the last time point as the final outcome of the
%time course - a bit more accurate
diff = [];
% NewRealradialAvgNuc.nucAvg(1,1)
dif = abs((NewRealradialAvgNuc{4} - SimradialAvgNuc{4}));
diff(:,1) = dif(:,3).^2;
costF = sum(sum(diff));
% costF = mean(mean(diff));
%     
display(["cost is", num2str(costF)])
end
