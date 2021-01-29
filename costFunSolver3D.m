
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
ic = [32 2 1.0000e-04];
matA = ic(1).*edge;
matI2 = ic(2).*edge;
matI3 = ic(2).*edge;

%9 parameters
% B1 = param.B1
K11_1 = param.K11_1 = vec(1)
K11_2 = param.K11_2 = vec(2)
n11 = param.n11 = vec(3)
kd1 = param.kd1 = vec(4)
K21 = param.K21 = vec(5)
% K41_1 = param.K41_1 = vec(2)
% K41_2 = param.K41_2 = vec(2)
% n41 = param.n41 = vec(2)
% B2 = param.B2 = vec(2)
K12_1 = param.K12_1 = vec(6)
K12_2 = param.K12_2 = vec(7)
n12 = param.n12 = vec(8)
K32 = param.K32 = vec(9)
kd2 = param.kd2 = vec(10)
K13_1 = param.K13_1 = vec(11)
K13_2 = param.K13_2 = vec(12)
n13 = param.n13 = vec(13)
% K23 = param.K23 = vec(2)
kd3 = param.kd3 = vec(14)
% kdm = param.kdm = vec(2)



for i = 2:2:8
    
    f = @(t,y) parametersfun(y, param);
    
    tspan = [i-2;i];
    sol = ode45(f,tspan,[matA(:);matI(:)]);
    y = sol.y(:,end);
    n = length(y)/3;
A = y(1:n);
I2 = y((n+1):2*n);
I3 = y((2*n+1):end);

matA = reshape(A,[],sqrt(n));
matI2 = reshape(I2,[],sqrt(n)); %lefty1?
matI3 = reshape(I3,[],sqrt(n)); %lefty2?
    t = round(sol.x(end));
    imshow(matI2,[])
    pause(1)
    imwrite(matA./65,[saveInPath,'A1_',num2str(t),'h.tif']);
    imwrite(matI2./65,[saveInPath,'I2_',num2str(t),'h.tif']);
    imwrite(matI3./65,[saveInPath,'I3_',num2str(t),'h.tif']);

end

%% FROM THE SIMULATIONS
dapi = imread('/Volumes/storage/Eleana/modelling_gastruloids/NODALmodellling/outPutODE45_2020/dapi0.tif');
dapi = dapi ./ 500;
dapi = imfill(dapi,'holes');

[SimradialAvgNuc]  = SimColoniesAnalysis(dapi);
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
% figure;plot(NewRealradialAvgNuc{2}(:,1));
figure;plot(SimradialAvgNuc{1}(:,1));
pause (2)
SimradialAvgNuc{1}(:,1)
% NewRealradialAvgNuc{2}(:,1)
dif1 = abs((NewRealradialAvgNuc{2}(:,1) - SimradialAvgNuc{1}(:,1))); % the 2nd condition of real col is 48h, the 1st for sim is 48h
dif2 = abs((NewRealradialAvgNuc{4}(:,1) - SimradialAvgNuc{3}(:,1))); % the 4nd condition of real col is 24h, the 3st for sim is ~24h

diff = dif1(:,1).^2 +dif2(:,1).^2 ; %(WNT channel is the first)
costF = sum(sum(diff));
% costF = mean(mean(diff));
%     
display(["cost is", num2str(costF)])
close all
end
