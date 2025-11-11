% Parameters
phi = 0.58;
alpha = 5/4;
g = 9.81;
theta = 25;
mu_b1 = tand(20.16);
mu_b2 = tand(37.65);
I0 =  0.434;
d = 0.001;
k = 1.0;
L = 5000 * d;
Nx = 50;
delta_x = L / (Nx - 1);
x = linspace(0, L, Nx);
T = 10;
Nt = 100;
delta_t = T / (Nt - 1);
t = linspace(0, T, Nt);

n = Nx;
m = Nt;
h = 39 * d * ones(n, m); % Initial values of h
u = 5 * ones(n, m); % Initial values of u
h_bc = 39 * d;
u_bc = 5;

% Calculate residuals
F = calculate_residuals(h, u, n, m, delta_t, delta_x, alpha, g, phi, theta, mu_b1, mu_b2, I0, d, k, h_bc, u_bc);

% Calculate Jacobian matrix
J = calculate_jacobian(h, u, n, m, delta_t, delta_x, alpha, g, phi, theta, mu_b1, mu_b2, I0, d, k, h_bc, u_bc);

% Newton-Raphson update
delta_X = -J^(-1) *F;


for i = 1:n-1
    for j = 1:m
        h(i+1, j) = h(i, j) +delta_X((i - 1) * m + j);
        u(i+1, j) = u(i, j) + delta_X(n * m + (i - 1) * m + j);
    end
end

% Display updated h and u
disp('Updated h:');

disp('Updated u:');
disp(u);

% After the loop for updating h and u:

% Specify the directory to save the plots
output_dir = 'Time';

% Create the folder if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Loop through all time steps (from t = 1 to T)
for t_idx = 1:m % m is the number of time steps
    % Extract the appropriate column for the current time step
    figure;
    
    % Plot h vs x
    subplot(2, 1, 1);
    plot(x, h(:, t_idx)); 
    xlabel('x');
    ylabel('h');
    title(['h vs x at t = ' num2str(t(t_idx)) ' sec']);
    
    % Plot u vs x
    subplot(2, 1, 2);
    plot(x, u(:, t_idx)); 
    xlabel('x');
    ylabel('u');
    title(['u vs x at t = ' num2str(t(t_idx)) ' sec']);
    
    % Save the figure in the specified folder
    saveas(gcf, fullfile(output_dir, ['h_u_t' num2str(t(t_idx)) '.png'])); % Save as PNG
end
%Calculating Residual Matrix Function

    function F = calculate_residuals(h, u, n, delta_t, delta_x, alpha, g, phi, theta, mu_b1, mu_b2, I0, d, k, h_bc, u_bc)
    F = zeros(2*n*n, 1);
    g_phi_cos_theta = g * phi * cosd(theta);
    idx = 1;
    for i = 1:n
        for j = 1:n
            if i-1 == 0
                h_im1 = h_bc;
                u_im1 = u_bc;
            else
                h_im1 = h(i-1, j);
                u_im1 = u(i-1, j);
            end
            if i+1 == n+1
                h_ip1 = h(i, j);
                u_ip1 = u(i, j);
            else
                h_ip1 = h(i+1, j);
                u_ip1 = u(i+1, j);
            end
            if j-1 == 0
                h_jm1 = h_bc;
                u_jm1 = u_bc;
            else
                h_jm1 = h(i, j-1);
                u_jm1 = u(i, j-1);
            end
            if j+1 == n+1
                h_jp1 = h(i, j);
                u_jp1 = u(i, j);
            else
                h_jp1 = h(i, j+1);
                u_jp1 = u(i, j+1);
            end
            
            hi = h(i,j);
            ui = u(i,j);
            Ib = (5 * d * ui) / (2 * hi * sqrt(g_phi_cos_theta * hi));
            mu_b = mu_b1 + (mu_b2 - mu_b1) /  (I0 /Ib)+1;
            
            F(idx) = h_jp1 - h_jm1 + 2*delta_t*(hi*(u_ip1 - u_im1)/(2*delta_x)) + 2*delta_t*(ui*(h_ip1 - h_im1)/(2*delta_x));
            idx = idx + 1;
            
            F(idx) = hi*(u_jp1 - u_jm1)/(2*delta_t) + ui*(h_jp1 - h_jm1)/(2*delta_t) ...
                     + alpha*(2*hi*ui*(u_ip1 - u_im1)/(2*delta_x) + ui^2*(h_ip1 - h_im1)/(2*delta_x)) ...
                     - g*hi*cosd(theta)*(tand(theta) - mu_b - k*(h_ip1 - h_im1)/(2*delta_x));
            idx = idx + 1;
        end
    end
end

%Calculating Jacobian Matrix Function

    function J = calculate_jacobian(h, u, n, delta_t, delta_x, alpha, g, phi, theta, mu_b1, mu_b2, I0, d, k, h_bc, u_bc)
    J = zeros(2*n*n, 2*n*n);
    g_phi_cos_theta = g * phi * cosd(theta);
    idx = 1;
    for i = 1:n
        for j = 1:n
            if i-1 == 0
                h_im1 = h_bc;
                u_im1 = u_bc;
            else
                h_im1 = h(i-1, j);
                u_im1 = u(i-1, j);
            end
            if i+1 == n+1
                h_ip1 = h(i, j);
                u_ip1 = u(i, j);
            else
                h_ip1 = h(i+1, j);
                u_ip1 = u(i+1, j);
            end
            if j-1 == 0
                h_jm1 = h_bc;
                u_jm1 = u_bc;
            else
                h_jm1 = h(i, j-1);
                u_jm1 = u(i, j-1);
            end
            if j+1 == n+1
                h_jp1 = h(i, j);
                u_jp1 = u(i, j);
            else
                h_jp1 = h(i, j+1);
                u_jp1 = u(i, j+1);
            end
            
            hi = h(i,j);
            ui = u(i,j);
            Ib = (5 * d * ui) / (2 * hi * sqrt(g_phi_cos_theta * hi));
            mu_b = mu_b1 + ((mu_b2 - mu_b1) / (I0/Ib))+1;
            
            % Calculate partial derivatives for F1
            dF1_dhij = 2*delta_t*(u_ip1 - u_im1)/(2*delta_x);
            dF1_dhijm1 = -1;
            dF1_dhijp1 = 1;
            dF1_dhim1j = -2*delta_t*(ui/(2*delta_x));
            dF1_dhip1j = 2*delta_t*(ui/(2*delta_x));
            dF1_duij = 2*delta_t*(h_ip1 - h_im1)/(2*delta_x);
            dF1_duim1j = -2*delta_t*(hi/(2*delta_x));
            dF1_duip1j = 2*delta_t*(hi/(2*delta_x));
            
            % Fill in Jacobian for F1
            J(idx, (i-1)*n + j) = dF1_dhij;
            if j-1 > 0
                J(idx, (i-1)*n + j-1) = dF1_dhijm1;
            end
            if j+1 <= n
                J(idx, (i-1)*n + j+1) = dF1_dhijp1;
            end
            if i-1 > 0
                J(idx, (i-2)*n + j) = dF1_dhim1j;
            end
            if i+1 <= n
                J(idx, i*n + j) = dF1_dhip1j;
            end
            J(idx, n*n + (i-1)*n + j) = dF1_duij;
            if i-1 > 0
                J(idx, n*n + (i-2)*n + j) = dF1_duim1j;
            end
            if i+1 <= n
                J(idx, n*n + i*n + j)= dF1_duip1j;
            end
            idx = idx + 1;
            
            % Calculate partial derivatives for F2
            dF2_dhij = (u_jp1 - u_jm1)/(2*delta_t) + alpha*(2*u(i,j)*(u_ip1 - u_im1)/(2*delta_x) + u(i,j)^2/(2*delta_x)) - g*cosd(theta)*(-mu_b - k*(h_ip1 - h_im1)/(2*delta_x));
            dF2_dhijm1 = -u(i,j)/(2*delta_t);
            dF2_dhijp1 = u(i,j)/(2*delta_t);
            dF2_dhim1j = -alpha*(2*u(i,j)*(u_ip1 - u_im1)/(2*delta_x) + u(i,j)^2/(2*delta_x));
            dF2_dhip1j = alpha*(2*u(i,j)*(u_ip1 - u_im1)/(2*delta_x) + u(i,j)^2/(2*delta_x));
            dF2_duij = hi*(h_jp1 - h_jm1)/(2*delta_t) + 2*alpha*(hi*(u_ip1 - u_im1)/(2*delta_x));
            dF2_duijm1 = -hi/(2*delta_t);
            dF2_duijp1 = hi/(2*delta_t);
            dF2_duim1j = -alpha*2*hi*u(i,j)/(2*delta_x);
            dF2_duip1j = alpha*2*hi*u(i,j)/(2*delta_x);
            
            % Fill in Jacobian for F2
            J(idx, (i-1)*n + j) = dF2_dhij;
            if j-1 > 0
                J(idx, (i-1)*n + j-1) = dF2_dhijm1;
            end
            if j+1 <= n
                J(idx, (i-1)*n + j+1) = dF2_dhijp1;
            end
            if i-1 > 0
                J(idx, (i-2)*n + j) = dF2_dhim1j;
            end
            if i+1 <= n
                J(idx, i*n + j) = dF2_dhip1j;
            end
            J(idx, n*n + (i-1)*n + j) = dF2_duij;
            if j-1 > 0
                J(idx, n*n + (i-1)*n + j-1) = dF2_duijm1;
            end
            if j+1 <= n
                J(idx, n*n + (i-1)*n + j+1) = dF2_duijp1;
            end
            if i-1 > 0
                J(idx, n*n + (i-2)*n + j) = dF2_duim1j;
            end
            if i+1 <= n
                J(idx, n*n + i*n + j) = dF2_duip1j;
            end
            idx = idx + 1;
        end
    end
end