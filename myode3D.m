function [dy] = myode3D(~,y) % a matrix 100X100
n = length(y)/3;
A = y(1:n);
I2 = y((n+1):2*n);
I3 = y((2*n+1):end);

matA = reshape(A,[],sqrt(n));
matI2 = reshape(I2,[],sqrt(n)); %lefty1?
matI3 = reshape(I3,[],sqrt(n)); %lefty2?
% %sapna's style

% dC =[2.0000e-04 2.0000e-04 0]
% B1 = 0.1000
K11_1 = 0.1100
K11_2 = 1
n11 = 5
kd1 = 1.0000e-03
K21 = 0
% K41_1 = 0
% K41_2 = 1
% n41 = 5
% B2 = 0.1000
K12_1 = 0.0200
K12_2 = 0.5000;
n12 = 5;
K32 = 0.0300;
kd2 = 1.0000e-03
K13_1 = 0.0100;
K13_2 = 12;
n13 = 5
% K23 = 0
kd3 = 2.0000e-04
% kdm = 1

% tmax = 1000
% dt = 1
% ic = [32 2 1.0000e-04]

% Da = 12; %or 10-12? 20 seems to work too
% ba = 0.2; %basal production rate of A
% bi = 0.1; %basal production rate of I
% ka  = 0.01; %inhibition of A on A
% kia = 0.02; %inhibition of I on A
% kda = 0.1; %degradation of A in the colony
% kdi = 0.1; %degradation of I in the colony
% kaa = 1.1;
% matAA = matA;
% matAA(matAA < 0.001) = 0 ;
% 
% matII = matI;
% matII(matII< 0.001) = 0 ;

% activator
matA = (K11_1).*(matA.^n11)./((1+(K21.*matI2)).*(1+K11_2.*matA.^n11))-(kd1).*matA;
%inhibitor
matI2 = (K12_1).*(matA.^n12)./(K12_2^n12+matA.^n12)-(kd2).*matI2-(K32).*matI3;%- matII.^3
 
matI3 = (K13_1).*matA.^n13./(K13_2.^n13+matA.^n13)-(kd3).*matI3;

dy = [matA(:);matI2(:);matI3(:)]
end