function [dy] = parametersfun(y, param)

%9 parameters
Da = param.Da; 
ba = param.ba; %basal production A 0.02
bi = param.bi; %basal production I 0.03
ka = param.ka; %autoproduction A 0.04
% kb = param.kb; % production due to BMP4 0.2
kia = param.kia; %inhibition due to I 0.05
kda = param.kda; % basal degradation A 0.015
kaa = param.kaa; %production of I 0.03
kdi = param.kdi; % basal degradation I 0.015
% B = 0.3; % concentration BMP4

n = length(y)/2;
A = y(1:n);
I = y((n+1):end);
matA = reshape(A,[],sqrt(n));
matI = reshape(I,[],sqrt(n));

% activator
matA = (ba+(matA)./((1+kia*matI).*(1+ka*matA)) - kda*matA + Da*del2(matA));
%inhibitor
matI = (kaa*matA + bi - kdi*matI );%- matII.^3 

% %new exp functions
% matA = (ba*B)./((kb+B).*(ka+matA)) - 1./(kia+matI) - kda*matA + Da.*del2(matA);
% matI = bi./(kaa+matA)- kdi*matI; %- matI.^3;

dy = [matA(:);matI(:)];

end