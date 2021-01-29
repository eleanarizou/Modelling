function [dy] = parametersfun3d(y, param)

%9 parameters
% dC =[2.0000e-04 2.0000e-04 0]
% B1 = param.B1
K11_1 = param.K11_1;
K11_2 = param.K11_2;
n11 = param.n11;
kd1 = param.kd1;
K21 = param.K21;
% K41_1 = param.K41_1
% K41_2 = param.K41_2
% n41 = param.n41
% B2 = param.B2
K12_1 = param.K12_1;
K12_2 = param.K12_2;
n12 = param.n12;
K32 = param.K32;
kd2 = param.kd2;
K13_1 = param.K13_1;
K13_2 = param.K13_2;
n13 = param.n13;
% K23 = param.K23
kd3 = param.kd3;
% kdm = param.kdm

n = length(y)/3;
A = y(1:n);
I2 = y((n+1):2*n);
I3 = y((2*n+1):end);

matA = reshape(A,[],sqrt(n));
matI2 = reshape(I2,[],sqrt(n)); %lefty1?
matI3 = reshape(I3,[],sqrt(n)); %lefty2?
% %sapna's style



% tmax = 1000
% dt = 1
% ic = [32 2 1.0000e-04]

% activator
matA = ((K11_1).*(matA.^n11)./((1+(K21.*matI2)).*(1+K11_2.*matA.^n11))-(kd1).*matA);
%inhibitor
matI2 = ((K12_1).*(matA.^n12)./(K12_2^n12+matA.^n12)-(kd2).*matI2-(K32).*matI3);%- matII.^3
 
matI3 = ((K13_1).*matA.^n13./(K13_2.^n13+matA.^n13)-(kd3).*matI3);

dy = [matA(:);matI2(:);matI3(:)];

end