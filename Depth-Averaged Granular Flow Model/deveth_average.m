clear all; close all; clc;
global h_initial vx_initial bin_h bin_t I_not mu_s mu_m d theta phi g alpha k h_initial0 vx_initial0 L_x
ts = 10; % Reduce time step

%% System variables
g = 9.81;
theta = 25;

%% System properties
d = 0.001;
L_x = 1000*d;

%% Rheological model parameters
I_not = 1;
mu_s = 0.1;
mu_m = 0.75;
phi = 0.58;
alpha = 5/4; % Bagnold velocity profile
k = 1; % Neglecting normal stress differences

%% Initial values
bin_h = 100;
bin_t = 100;
h_initial =39*d*ones(1, bin_h);
vx_initial = 5*ones(1, bin_h);
h_initial(1) = 39* d;
vx_initial(1) = 5;
h_initial0 = 39* d;
vx_initial0 = 5;

%% Initialization
y = linspace(0, L_x, bin_h);
t = linspace(0, ts, bin_t);

count = 1;
time = 0;

% ODE solver options
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 1e-2); % Adjusted tolerances and max step

output_dir = 'Time';
if ~exist(output_dir, 'dir')
    status = mkdir(output_dir);
    if status == 0
        error('Failed to create directory: %s', output_dir);
    else
        disp(['Directory created: ', output_dir]);
    end
else
    disp(['Directory already exists: ', output_dir]);
end

%% Simulation loop
while time < 500
    m = 0; % Symmetry parameter for PDEPE

    % Solve the PDE
    try
        sol = pdepe(m, @pdex2pde, @pdex2ic, @pdex2bc, y, t, options); 
    catch ME
        warning('Time integration has failed. Solution is available at requested time points up to t=%.2f.\nError: %s', time, ME.message);
        break;
    end

    % Extract solutions
    velocity = sol(:,:,1);
    height = sol(:,:,2);

    count = count + 1;
    time = time + ts;
    disp(['Current simulation time: ', num2str(time)]);

    %% Properties calculation
    vx = velocity(end, :); % Update velocity, vx after time step ts
    h = height(end, :);

    %% Restricting f to be in between 0 to 1
    h_initial = h;
    vx_initial = vx;

    % Update data
    data = [y', vx', h'];

    % Write data to file
    filename = fullfile(output_dir, sprintf('properties-%05d.txt', count)); % Adjusted filename format
    disp(['Writing data to file: ', filename]);
    try
        writematrix(data, filename, 'Delimiter', ' ');
    catch ME
        error('Failed to write file: %s\nError: %s', filename, ME.message);
    end

    % Plot and save the figure
    figure;
    subplot(2, 1, 1);
    plot(y, h,'linewidth',2);
    title(sprintf('Height at time = %.2f', time));
    xlabel('Position y');
    ylabel('Height h');

    subplot(2, 1, 2);
    plot(y, vx,'linewidth',2);
    title(sprintf('Velocity at time = %.2f', time));
    xlabel('Position y');
    ylabel('Velocity vx');

    saveas(gcf, fullfile(output_dir, sprintf('plot-%05d.png', count)));
    close
end
function [pl ,ql ,pr ,qr] = pdex2bc(xl ,ul ,xr ,ur ,t)
    global time h_initial0 vx_initial0
    pl = [ul(1) - vx_initial0; ul(2) - h_initial0];
    % pl = [ul(1); ul(2)];
    ql = [0; 0];

    pr = [0; 0];
    qr = [1; 1]; 
end
function u0 = pdex2ic(x)
    global bin_h L_x vx_initial h_initial time
  
    s = (x - 0) / ((L_x) / (bin_h - 1)) + 1;
    s = round(s); % bin index
  
    f0 = h_initial(s);
    v0 = vx_initial(s);
    u0 = [v0; f0];
end
function [C, F, S] = pdex2pde(x, t, u, DuDx)
    global I_not mu_s mu_m d theta phi g alpha k
  
    vx = u(1);
    h = u(2);
    dvx_dy = DuDx(1);
    dh_dy = DuDx(2);

    I = 5 * d * vx / (2 * h * sqrt(phi * g * h * cosd(theta)));
    
    mu_b = mu_s + (mu_m - mu_s) / (1 + I_not / I); 
   
    C1 = h ;
    S1 =  g * h * cosd(theta) * (tand(theta) - mu_b) - k * g * h * cosd(theta) * dh_dy+vx^2*(1-alpha)*dh_dy+vx*h*(1-2*alpha)*dvx_dy;
    F1 = 0;

    C2 = 1;
    S2 = -h*dvx_dy-vx*dh_dy;
    F2 = 0;

    C = [C1; C2];
    S = [S1; S2];
    F = [F1; F2];
end
