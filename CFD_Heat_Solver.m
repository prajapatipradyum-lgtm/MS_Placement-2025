%% Heat Equation Solver

clear all; close all; clc;

fprintf('=== DEBUGGED HEAT EQUATION SOLVER ===\n');

%% ============================================================================
%                     1D HEAT EQUATION 
%% ============================================================================

% Parameters
L = 1.0;                    % Length of rod (m)
alpha = 1e-3;               % Thermal diffusivity (m²/s)
T_left = 100;               % Left boundary temperature (°C)
T_right = 0;                % Right boundary temperature (°C)
T_initial = 20;             % Initial temperature (°C)

% Numerical Parameters
nx = 51;                    % Number of spatial grid points
dx = L / (nx - 1);          % Spatial step size
dt = 0.1;                   % Time step

% Simulation time
t_char = L^2 / alpha;       % Characteristic time
t_final = 2 * t_char;       % Simulation time
nt = round(t_final / dt);   

% Stability check
CFL = alpha * dt / dx^2;
fprintf('CFL number: %.4f\n', CFL);

% Initialize arrays
x = linspace(0, L, nx);
T = ones(nx, 1) * T_initial;

% Apply initial boundary conditions
T(1) = T_left;
T(end) = T_right;

% DEBUG: Print initial conditions
fprintf('Initial conditions:\n');
fprintf('  T(1) = %.1f°C, T(end) = %.1f°C, T(center) = %.1f°C\n', ...
        T(1), T(end), T(round(nx/2)));

% Storage for snapshots
snapshot_times = [0, t_char/4, t_char/2, t_char, t_final];
snapshot_indices = round(snapshot_times / dt) + 1;
T_snapshots = zeros(nx, length(snapshot_times));
T_snapshots(:, 1) = T;  % Store initial condition

fprintf('Starting 1D simulation: %d time steps\n', nt);

%% Time stepping loop
snapshot_counter = 2;
T_old = T;  % Initialize T_old

for n = 1:nt
    % Store current state before updating
    T_old(:) = T;
    
    % Update interior points only
    for i = 2:nx-1
        T(i) = T_old(i) + alpha * dt / dx^2 * ...
               (T_old(i+1) - 2*T_old(i) + T_old(i-1));
    end
    
    % CRITICAL: Enforce boundary conditions after each update
    T(1) = T_left;
    T(end) = T_right;
    
    % DEBUG: Check temperature range periodically
    if mod(n, 1000) == 0
        fprintf('  Step %d: T_range=[%.1f,%.1f], T_center=%.1f°C\n', ...
                n, min(T), max(T), T(round(nx/2)));
    end
    
    % Store snapshots
    if snapshot_counter <= length(snapshot_indices) && n == snapshot_indices(snapshot_counter)
        T_snapshots(:, snapshot_counter) = T;
        fprintf('Snapshot %d at t=%.1fs: T_center=%.1f°C\n', ...
                snapshot_counter, (n-1)*dt, T(round(nx/2)));
        snapshot_counter = snapshot_counter + 1;
    end
end

% DEBUG: Final temperature check
fprintf('Final 1D state after time loop:\n');
fprintf('  T(1)=%.1f, T(end)=%.1f, T(center)=%.1f°C\n', ...
        T(1), T(end), T(round(nx/2)));
fprintf('  min(T)=%.1f, max(T)=%.1f\n', min(T), max(T));

%% Store final snapshot (CRITICAL FIX)
if snapshot_counter <= length(snapshot_indices)
    T_snapshots(:, end) = T;  % Store final state
    fprintf('Final snapshot stored: T_center=%.1f°C\n', T_snapshots(round(nx/2), end));
end

%% ============================================================================
%                     2D HEAT EQUATION (SIMPLIFIED)
%% ============================================================================

fprintf('\n=== 2D HEAT EQUATION ===\n');

% 2D Parameters
Lx = 1.0; Ly = 1.0;
alpha_2d = 1e-3;
T_top = 100; T_bottom = 0; T_left_2d = 50; T_right_2d = 50;
T_initial_2d = 25;

% Grid
nx_2d = 21; ny_2d = 21;
dx_2d = Lx / (nx_2d - 1);
dy_2d = Ly / (ny_2d - 1);
dt_2d = 0.01;

% Shorter simulation time for testing
t_final_2d = 500;  % Reduced time
nt_2d = round(t_final_2d / dt_2d);

% Initialize
x_2d = linspace(0, Lx, nx_2d);
y_2d = linspace(0, Ly, ny_2d);
[X, Y] = meshgrid(x_2d, y_2d);

T_2d = ones(ny_2d, nx_2d) * T_initial_2d;
T_2d_old = T_2d;

% Apply boundary conditions
T_2d(1, :) = T_bottom;
T_2d(end, :) = T_top;
T_2d(:, 1) = T_left_2d;
T_2d(:, end) = T_right_2d;

fprintf('2D simulation: %d time steps\n', nt_2d);

%% 2D Time stepping
for n = 1:nt_2d
    T_2d_old(:,:) = T_2d;
    
    % Update interior points
    for i = 2:ny_2d-1
        for j = 2:nx_2d-1
            T_2d(i,j) = T_2d_old(i,j) + alpha_2d * dt_2d * ...
                       ((T_2d_old(i+1,j) - 2*T_2d_old(i,j) + T_2d_old(i-1,j)) / dy_2d^2 + ...
                        (T_2d_old(i,j+1) - 2*T_2d_old(i,j) + T_2d_old(i,j-1)) / dx_2d^2);
        end
    end
    
    % Enforce 2D boundary conditions
    T_2d(1, :) = T_bottom;
    T_2d(end, :) = T_top;
    T_2d(:, 1) = T_left_2d;
    T_2d(:, end) = T_right_2d;
    
    % Progress
    if mod(n, round(nt_2d/5)) == 0
        center_temp = T_2d(round(ny_2d/2), round(nx_2d/2));
        fprintf('2D Progress %.0f%%: Center T = %.1f°C\n', 100*n/nt_2d, center_temp);
    end
end

%% ============================================================================
%                          ANALYSIS WITH DEBUGGING
%% ============================================================================

% Analytical solution
T_analytical = T_left - (T_left - T_right) * x / L;

% DEBUG: Check analytical solution
fprintf('\nAnalytical solution check:\n');
fprintf('  T_analytical(1) = %.1f°C\n', T_analytical(1));
fprintf('  T_analytical(center) = %.1f°C\n', T_analytical(round(nx/2)));
fprintf('  T_analytical(end) = %.1f°C\n', T_analytical(end));

% DEBUG: Check final snapshot
fprintf('\nFinal snapshot check:\n');
fprintf('  T_snapshots size: [%d x %d]\n', size(T_snapshots));
fprintf('  T_snapshots(1, end) = %.1f°C\n', T_snapshots(1, end));
fprintf('  T_snapshots(center, end) = %.1f°C\n', T_snapshots(round(nx/2), end));
fprintf('  T_snapshots(end, end) = %.1f°C\n', T_snapshots(end, end));

% Calculate error properly
center_idx = round(nx/2);
T_final_numerical = T_snapshots(:, end);  % Final numerical solution
error_1d = abs(T_final_numerical - T_analytical');
max_error_1d = max(error_1d);

% Get individual values for validation
T_center_numerical = T_snapshots(center_idx, end);
T_center_analytical = T_analytical(center_idx);
T_min_final = min(T_final_numerical);
T_max_final = max(T_final_numerical);

% 2D results
center_T_2d = T_2d(round(ny_2d/2), round(nx_2d/2));

%% PLOTTING
figure(1);
set(gcf, 'Position', [100 100 1200 800]);

% 1D Evolution plot
subplot(2,3,1);
colors = lines(length(snapshot_times));
for i = 1:length(snapshot_times)
    plot(x, T_snapshots(:,i), 'Color', colors(i,:), 'LineWidth', 2);
    hold on;
end
plot(x, T_analytical, 'k--', 'LineWidth', 3);
xlabel('Position x (m)');
ylabel('Temperature (°C)');
title('1D Temperature Evolution');
legend([cellstr(num2str(snapshot_times', 't=%.0fs')); {'Analytical'}], 'Location', 'best');
grid on;

% Final comparison
subplot(2,3,2);
plot(x, T_final_numerical, 'r-', 'LineWidth', 2);
hold on;
plot(x, T_analytical, 'k--', 'LineWidth', 2);
xlabel('Position x (m)');
ylabel('Temperature (°C)');
title('1D Final vs Analytical');
legend('Numerical', 'Analytical');
grid on;

% Error plot
subplot(2,3,3);
plot(x, error_1d, 'b-o', 'MarkerSize', 4);
xlabel('Position x (m)');
ylabel('Absolute Error (°C)');
title(sprintf('1D Error (max=%.3f°C)', max_error_1d));
grid on;

% 2D results
subplot(2,3,4);
contourf(X, Y, T_2d, 15);
colorbar;
title('2D Temperature Field');
xlabel('x (m)'); ylabel('y (m)');

subplot(2,3,5);
surf(X, Y, T_2d);
shading interp;
colorbar;
title('2D Temperature Surface');
xlabel('x (m)'); ylabel('y (m)'); zlabel('T (°C)');

% Cross-sections
subplot(2,3,6);
center_row = round(ny_2d/2);
center_col = round(nx_2d/2);
plot(x_2d, T_2d(center_row, :), 'r-', 'LineWidth', 2);
hold on;
plot(y_2d, T_2d(:, center_col), 'b-', 'LineWidth', 2);
xlabel('Position (m)');
ylabel('Temperature (°C)');
title('2D Cross-sections');
legend('Horizontal', 'Vertical');
grid on;

%% FINAL VALIDATION 
fprintf('\n=== DEBUGGED VALIDATION RESULTS ===\n');

fprintf('1D Heat Equation:\n');
fprintf('  • Max error vs analytical: %.3f °C\n', max_error_1d);
fprintf('  • Final temp range: [%.1f, %.1f] °C\n', T_min_final, T_max_final);
fprintf('  • Center temp: %.1f °C (analytical: %.1f °C)\n', T_center_numerical, T_center_analytical);
fprintf('  • Relative error at center: %.2f%%\n', 100*abs(T_center_numerical-T_center_analytical)/T_center_analytical);

fprintf('2D Heat Equation:\n');
fprintf('  • Center temperature: %.1f °C (change: %.1f°C)\n', center_T_2d, center_T_2d - T_initial_2d);
fprintf('  • Temp range: [%.1f, %.1f] °C\n', min(T_2d(:)), max(T_2d(:)));

% Validation checks
fprintf('\nValidation Summary:\n');
success_count = 0;

if max_error_1d < 5.0
    fprintf('✅ 1D numerical error acceptable (%.3f°C < 5°C)\n', max_error_1d);
    success_count = success_count + 1;
else
    fprintf('❌ 1D error too high: %.1f°C\n', max_error_1d);
end

if abs(T_center_numerical - T_center_analytical) < 2.0
    fprintf('✅ 1D center temperature accurate (error: %.1f°C)\n', abs(T_center_numerical - T_center_analytical));
    success_count = success_count + 1;
else
    fprintf('❌ 1D center temperature error: %.1f°C\n', abs(T_center_numerical - T_center_analytical));
end

if abs(T_min_final - T_right) < 0.1 && abs(T_max_final - T_left) < 0.1
    fprintf('✅ 1D boundary conditions maintained\n');
    success_count = success_count + 1;
else
    fprintf('❌ 1D boundary conditions: min=%.1f (should be %.1f), max=%.1f (should be %.1f)\n', ...
            T_min_final, T_right, T_max_final, T_left);
end

if abs(center_T_2d - T_initial_2d) > 5
    fprintf('✅ 2D shows significant heat transfer (change: %.1f°C)\n', abs(center_T_2d - T_initial_2d));
    success_count = success_count + 1;
else
    fprintf('❌ 2D minimal heat transfer (change: %.1f°C)\n', abs(center_T_2d - T_initial_2d));
end

fprintf('\nSUCCESS RATE: %d/4 tests passed\n', success_count);
if success_count >= 3
    fprintf('🎉 CFD PROJECT SUCCESSFULLY COMPLETED! 🎉\n');
else
    fprintf('⚠️  Project needs further debugging\n');
    
    % Additional debugging info
    fprintf('\nDEBUG INFO:\n');
    fprintf('  • Final T array size: %d\n', length(T));
    fprintf('  • Snapshots array size: [%d x %d]\n', size(T_snapshots));
    fprintf('  • Number of time steps completed: %d\n', nt);
    fprintf('  • Current T values: min=%.1f, max=%.1f, center=%.1f\n', min(T), max(T), T(center_idx));
end