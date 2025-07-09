% Parameters
clear
tic
D_D = 4;        % Diffusion coefficient for D
D_W = 3;        % Diffusion coefficient for W
D_S = 0;       % Diffusion coefficient for S
D_O = 0;       % Diffusion coefficient for O

s_WD = 1.1;       % Source term coefficient
s_W = .1;       % Source term coefficient
s_WS = 1.1;
a_W = 1;        % Parameter in the source term
a_D = 0.2;
gamma_DW = 0.1;   % Decay term
gamma_SO = 0.5;
gamma_WO = 0.5;
gamma_DS = 0.1;
k_D = 0.1;        % Decay coefficient for D
k_W = 0.1;        % Decay coefficient for W
k_S = 0.1;
k_O = 0.1;
k = 4;
K_D = 0.1;        % Saturation constant for D
K_W = 0.1;
K_S = 0.1;
sigma_W = .01;    % Source rate for D
sigma_O = .01;
W_threshold = 2;  % Threshold WNT concentration where SOX9 starts decreasing


% Spatial domain
L = 1000;          % Length of the spatial domain
x = linspace(0, L, 1000);  % Discretize the spatial domain

% Time domain
t_max = 4000;
t = linspace(0, t_max, 4000); % Time points for the solution

% Solve the PDE system using pdepe
m = 0;  % Symmetry parameter (0 for slab geometry in 1D)
sol = pdepe(m, @(x,t,u,DuDx) pde_system(x,t,u,DuDx,D_D,D_W,D_S,D_O,s_W,s_WD,s_WS,a_W,K_D,gamma_DW,gamma_SO,gamma_WO,gamma_DS,k_D,k_W,k_S,sigma_W,sigma_O,K_W,k_O,K_S, k), @initial_conditions, @boundary_conditions, x, t);
% Extract the solutions for D and W from the pdepe output
D_solution = sol(:,:,1); % Solution for D
W_solution = sol(:,:,2); % Solution for W
S_solution = sol(:,:,3); % Solution for D
O_solution = sol(:,:,4); % Solution for W

% Find the minimum and maximum values across both D_solution and W_solution
min_value = min([D_solution(:); W_solution(:); S_solution(:); O_solution(:)]);
max_value = max([D_solution(:); W_solution(:); S_solution(:); O_solution(:)]);

% Plot the solution for D over time
figure;
subplot(2,2,1);
surf(x, t, D_solution, 'EdgeColor', 'none');
title('Solution for D(x,t)');
xlabel('x');
ylabel('t');
zlabel('Dkk1');
view(2);  % View as a 2D plot with color mapping
caxis([min_value, max_value]);  % Set the color axis limits
colorbar;  % Add a colorbar to visualize the color scale
axis square
% Plot the solution for W over time
subplot(2,2,2);
surf(x, t, W_solution, 'EdgeColor', 'none');
title('Solution for Wnt');
xlabel('x');
ylabel('t');
zlabel('Wnt');
view(2);  % View as a 2D plot with color mapping
caxis([min_value, max_value]);  % Set the color axis limits to match D
colorbar;  % Add a colorbar to visualize the color scale
axis square
subplot(2,2,3);
surf(x, t, S_solution, 'EdgeColor', 'none');
title('Solution for SOX9');
xlabel('x');
ylabel('t');
zlabel('SOX9');
view(2);  % View as a 2D plot with color mapping
caxis([min_value, max_value]);  % Set the color axis limits
colorbar;  % Add a colorbar to visualize the color scale
axis square
% Plot the solution for W over time
subplot(2,2,4);
surf(x, t, O_solution, 'EdgeColor', 'none');
title('Solution for OTX2');
xlabel('x');
ylabel('t');
zlabel('OTX2');
view(2);  % View as a 2D plot with color mapping
caxis([min_value, max_value]);  % Set the color axis limits to match D
colorbar;  % Add a colorbar to visualize the color scale
axis square
colormap jet;

% time propagation movie %
create_wave_movie(W_solution, D_solution, S_solution, O_solution, t, x);


toc

function [c,f,s] = pde_system(x,t,u,DuDx,D_D,D_W,D_S,D_O,s_W,s_WD,s_WS,a_W,K_D,gamma_DW,gamma_SO,gamma_WO,gamma_DS,k_D,k_W,k_S,sigma_W,sigma_O,K_W,k_O,K_S,k)
    D = u(1); % DKK1 variable
    W = u(2); % WNT variable
    S = u(3); % SOX9 variable
    O = u(4); % OTX2 variable

    % Define thresholds for WNT accumulation and SOX9 formation
    W_threshold = 1.5;  % WNT threshold for SOX9 behavior change
    W_accumulation_threshold = 2.5;  % WNT threshold for accumulation

    % Modify WNT's diffusion coefficient based on DKK1 concentration
    if D > 1.0  % Threshold for DKK1 to start limiting WNT diffusion
        D_W_eff = D_W / (1 + D);  % Reduce WNT diffusion as DKK1 concentration increases
    else
        D_W_eff = D_W;  % Normal diffusion coefficient for WNT if DKK1 is low
    end

    % PDE for D (DKK1)
    c(1) = 1;  % Time derivative term coefficient for D
    f(1) = D_D * DuDx(1);  % Spatial derivative term for D (faster diffusion)
    s(1) = s_WD * W^2 / (1 + a_W * W) - gamma_DW * D - k_D * D;  % Source term for D (DKK1)

    % PDE for W (WNT), diffusion limited by DKK1 (D)
    c(2) = 1;  % Time derivative term coefficient for W
    f(2) = D_W_eff * DuDx(2);  % Spatial derivative term for W with reduced diffusion if D is high
    
    % WNT is inhibited by DKK1 (D) and SOX9 (S)
    s(2) = (s_W * W^2 / ((1 + a_W * W^2) + (K_D + D^2))) / (1 + D^2) - k_W * W * W + sigma_W;  % Source term for W (WNT)

    % PDE for SOX9 (S)
    c(3) = 1;  % Time derivative term coefficient for S
    f(3) = D_S * DuDx(3);  % Spatial derivative term for S
    s(3) = s_WS * W^2 / ((1 + a_W * W^2) + (K_S + S)) - gamma_SO * S * O - gamma_DS * D * S - k_S * S;
        
    % PDE for OTX2 (O), with SOX9 inhibiting OTX2
    c(4) = 1;  % Time derivative term coefficient for O
    f(4) = D_O * DuDx(4);  % Spatial derivative term for O
    
    % OTX2 is inhibited by SOX9, and WNT influences its behavior
    s(4) = -k_O * O + sigma_O / (K_W + W^2) - gamma_WO * W * O - k * O;

    % Ensure c, f, and s are column vectors
    c = [c(1); c(2); c(3); c(4)];
    f = [f(1); f(2); f(3); f(4)];
    s = [s(1); s(2); s(3); s(4)];
end



% Nested function to define the initial conditions
function u0 = initial_conditions(x)
    % Initialize D(x,0) = 0 for all x
    D0 = 0;

    % Initialize W(x,0) = non-zero near the edges (for example, at the first and last few cells)
    W0 = 0;  % Default value for W
    if x <= 1  % Set W to non-zero only at the first and last few edge cells
        W0 = 1;  % Example value for W at the edges
    end

    % Initialize D(x,0) = 0 for all x
    S0 = 1;

    % Initialize D(x,0) = 0 for all x
    O0 = 0;

    % Combine initial conditions for D and W
    u0 = [D0; W0; S0; O0];  % Return initial conditions as a column vector
end

% Nested function to define the boundary conditions (Neumann: no flux)
    function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
        % Dirichlet boundary conditions: D(0,t) = D(L,t) = 0, W(0,t) = W(L,t) = 0
        pl = [0; 0; 0; 0];  % Left boundary
        ql = [1; 1; 1; 1];  % Coefficient for Dirichlet BCs
        pr = [0; 0; 0; 0];  % Right boundary
        qr = [1; 1; 1; 1];  % Coefficient for Dirichlet BCs
    end

% % Nested function to define the boundary conditions (Dirichlet: fixed conc)
%     function [pl,ql,pr,qr] = boundary_conditions(xl,ul,xr,ur,t)
%         % Dirichlet boundary conditions: D(0,t) = D(L,t) = 0, W(0,t) = W(L,t) = 0
%         pl = [ul(1)-0; ul(2)-0; ul(3)-0; ul(4)-0];  % Left boundary
%         ql = [0; 0; 0; 0];  % Coefficient for Dirichlet BCs
%         pr = [ur(1)-0; ur(2)-0; ur(3)-0; ur(4)-0];  % Right boundary
%         qr = [0; 0; 0; 0];  % Coefficient for Dirichlet BCs
%     end

    function create_wave_movie(W_history, D_history, S_history, O_history, tspan, x)
    % Create figure for time-lapse movie
    figure;
    plot_W = plot(x, W_history(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Wnt');
    hold on;
    plot_D = plot(x, D_history(1,:), 'r--', 'LineWidth', 2, 'DisplayName', 'Dkk1');
    plot_S = plot(x, S_history(1,:), 'k--', 'LineWidth', 2, 'DisplayName', 'SOX9');
    plot_O = plot(x, O_history(1,:), 'g--', 'LineWidth', 2, 'DisplayName', 'OTX2');
   
    legend('show');
    xlabel('Position');
    ylabel('Concentration');
    title('Wnt, Dkk1, SOX9, OTX2 Propagation over Time');
%     ylim([0,10]);
    grid on;

    % Movie creation
    video_file = VideoWriter('/Users/elenirea/Documents/WarmflashWNTproject/WNT_turingfunction/wntmodelatrayee/Wnt_Dkk1_SOX9_OTX2_propagation.mp4', 'MPEG-4');
    open(video_file);

    for t_idx = 1:length(tspan)
        set(plot_W, 'YData', W_history(t_idx,:));  % Update Wnt plot
        set(plot_D, 'YData', D_history(t_idx,:));  % Update  Dkk1 plot
        set(plot_S, 'YData', S_history(t_idx,:));  % Update  SOX9 plot
        set(plot_O, 'YData', O_history(t_idx,:));  % Update  OTX2 plot
        drawnow;

        % Capture frame for movie
        frame = getframe(gcf);
        writeVideo(video_file, frame);
    end

    % Close video
    close(video_file);
end

