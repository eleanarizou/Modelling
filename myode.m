

function [dy] = myode(~,y) % a matrix 100X100
n = length(y)/2;
A = y(1:n);
I = y((n+1):end);
matA = reshape(A,[],sqrt(n));
matI = reshape(I,[],sqrt(n));

% %sapna's style

Da = 12; %or 10-12? 20 seems to work too
ba = 0.2; %basal production rate of A
bi = 0.1; %basal production rate of I
ka  = 0.01; %inhibition of A on A
kia = 0.02; %inhibition of I on A
kda = 0.1; %degradation of A in the colony
kdi = 0.1; %degradation of I in the colony
kaa = 1.1;
% matAA = matA;
% matAA(matAA < 0.001) = 0 ;
% 
% matII = matI;
% matII(matII< 0.001) = 0 ;

% activator
matA = (ba+(matA)./((1+kia*matI).*(1+ka*matA)) - kda*matA + Da*del2(matA));
%inhibitor
matI = (kaa*matA + bi - kdi*matI );%- matII.^3
   
dy = [matA(:);matI(:)]


%%my simpler style
% k1 = 0.05;
% k2 = 0.01;
% k3 = 0.04;
% k4 = 0.04;
% k5 = 0.06
% Pw = 5;
% Pa = 3;
% Dw = 10;
% 
% %WNT 
% 
% % ba + (matA.^2)./((1+k5/k1*matI).*(1+k4/k3*matA.^2))
% % matA = tau(1).*(k3*matA.^2 - k2*matI -k1*matA - k4*matA + Dw*del2(matA));  % - sum of second order partial derivatives in space  
% % matA = tau(1).*(Pw + (matA.^2)./((1+(k5/k1)*matI).*(1+(k4/k3)*matA.^2)) - k4*matA + Dw*del2(matA));  % - sum of second order partial derivatives in space  
% 
% matA = tau(1).*(Pw + k3*matA - k2*matI -k1*matA - k4*matA + Dw*del2(matA));  % - sum of second order partial derivatives in space  
% 
% %AXIN2 
% matI = tau(2).*(Pa + k1*matA.^2 - k5*matI ); %no diffusion- matII



end

