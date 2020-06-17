

saveInPath = '/Users/elenirea/Documents/WarmflashWNTproject/WNT_turingfunction/outPutODE45_2020/'
radius_colony = 600; %%Radius of colony
dtv = [10 120]; %time step maybe a smaller
dt = 5;
L= radius_colony*1.7067;  %domain size
N = 100;
x =linspace(-L,L,N);% (L/N)*(1:N)'; %find the coordinates for an image size of 512
y = x;
[X,Y] = ndgrid(x,x);

%in this area of code I can also define limits of certain patterns

%%%Define cutoff for colony
R = sqrt(X.^2 + Y.^2);
chi = ((1-tanh((R - radius_colony)/0.1))/2); %whatever positive plot , trash the rest , thresholding

%activator
edgewidth = 100;
edge = (1+tanh((R - radius_colony + edgewidth)/0.1))/2.*(1-tanh((R - radius_colony )/0.1));%.*(1+tanh((R - rad_col)/0.1));
media = ~chi;
%inhibitor
edgewidth = 200;
ed = (1+tanh((R - radius_colony + edgewidth)/0.1))/2.*(1-tanh((R - radius_colony )/0.1));%.*(1+tanh((R - rad_col)/0.1));


% i559 = imread('/Users/elenirea/Desktop/Image0448_01_559.tif.png')> 10000;
% i559 = imresize(i559, [100,100]);
% i635 = imread('/Users/elenirea/Desktop/Image0448_01_635.tif.png')> 10000;
% i635 = imresize(i635, [100,100]);
% im = [i559,i635];
% imshow(im)

% k = [0:N/2-1, 0 , -N/2+1:-1]'/(L/(2*pi));  %%fourier vector, for matlab fft - sets up a line that repeat for certain length and then reverse and repaet again for the same length-for mesh grid
% %set up the mesh grid?
% [xi,eta] = ndgrid(k,k);
% d = [6;0]
% %spdiags(B,d,m,n) creates an m-by-n sparse matrix from the columns of B and
% %places them along the diagonals specified by d.
% %it is not used?
% % D = spdiags(d',0,2,2); %%Sparse matrix formed from diagonals - the non zeros are compressed -faster
% %construct matrices L
% L1 = -d(1)*(eta.^2+xi.^2); % length for first morphogen?
% L2 = -d(2)*(eta.^2+xi.^2);
%
% % figure;
% % imagesc(L1);
% % figure;
% % imagesc(L2);
%
% %are these the boundary conditions?
% S1 = 1 - dt*L1;
% S2 = 1 - dt*L2;



Da = 12; %or 10-12? 20 seems to work too
ba = 0.2; %basal production rate of A
bi = 0.1; %basal production rate of I
ka  = 0.01; %inhibition of A on A
kia = 0.02; %inhibition of I on A
kda = 0.1; %degradation of A in the colony
kdi = 0.1; %degradation of I in the colony
kaa = 1.1;

f = @(t,y)parametersfun(y, Da, ba, bi, ka, kia, kda, kdi, kaa);

tic
tspan = [0;11];
ic = [0.5, 0.2];
sol = ode45(f,tspan,[ic(1).*edge(:);ic(2).*edge(:)]);
toc

% f1 = @(A,I) (ba +(A))./((1+kia*I).*(1+ka*A)) - kda*A +Da*del2(A);
% f2 = @(A,I) (kaa*A - kdi*I + bi);
%
% fff =@(a) [f1(a(1),a(2)) ;f2(a(1),a(2))];
% [equil] = fsolve(fff,[ic(1) ,ic(2)]);
% equil;

FigH = figure('Position', get(0, 'Screensize'));
F    = getframe(FigH);
for i = 1:length(sol.x)
    t = sol.x(i);
    y = sol.y(:,i);
    n = length(y)/2;
    A = y(1:n);
    I = y((n+1):end);
    matA = reshape(A,[],sqrt(n));
    matI = reshape(I,[],sqrt(n));
    mat = [matA,matI];
    if (i == 28 ) || (i == 36) ||(i == 49) ||(i == 63) ||(i == 72)
        subtightplot(1,7,round(i/7)-3,0.05)
        imagesc(2*x,2*x,mat);
        title(['A, I', num2str(sol.x(i))]);
        colorbar
        pause(1)
        %         t = round(sol.x(i));
        %         imwrite(matA./65,[saveInPath,'u1_',num2str(t),'h.tif']);
        %         imwrite(matI./65,[saveInPath,'u2_',num2str(t),'h.tif']);
    end
    %     close all
end


matA = ic(1).*edge
matI = ic(2).*edge
for i = 2:2:8

    f = @(t,y)parametersfun(y, Da, ba, bi, ka, kia, kda, kdi, kaa);

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
    imwrite(matA./65,[saveInPath,'u1_',num2str(t),'h.tif']);
    imwrite(matI./65,[saveInPath,'u2_',num2str(t),'h.tif']);
end



% for i = 1:length(sol.x)
%     t = sol.x(i);
%     y = sol.y(:,i);
%     n = length(y)/2;
%     A = y(1:n);
%     I = y((n+1):end);
%     matA = reshape(A,[],sqrt(n));
%     matI = reshape(I,[],sqrt(n));
%     if (i == 28 ) || (i == 36) ||(i == 49) ||(i == 63) ||(i == 72)
%         t = round(sol.x(i));
%         imwrite(matA./65,[saveInPath,'u1_',num2str(t),'h.tif']);
%         imwrite(matI./65,[saveInPath,'u2_',num2str(t),'h.tif']);
%     end
%     mat = [matA,matI];
%     figure;
%     imagesc(2*x,2*x,mat);
%     title(['A, I', num2str(sol.x(i))]);
%     colorbar
%     pause(1)
%     close all
% end


%% Save the whole sensitivity analysis - parameter values and result and plots

%     selected = split(ls(saveInPath));
%     pattern = ["20h","25h", "30h","40h","50h"];
%     timepoints = selected(contains(selected,pattern));
%

% Version='1';
% FileNameAndLocation=[mfilename('fullpath')];
% newbackup=sprintf('%sbackup%s.m',FileNameAndLocation,Version);
% currentfile=strcat(FileNameAndLocation, '.m');
% copyfile(currentfile,newbackup);
%
% save(['result_data_run_number_' num2str(n) ',mat'],'s')