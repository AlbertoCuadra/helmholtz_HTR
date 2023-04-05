% Routine to compute the dissipation, dissipation rate, and the solenoidal
% and compressive parts of a three-dimensional velocity field of a DNS
% obtained using the Hypersonic Task-based Research (HTR) solver [1]
% (see the GitHub repository [2])
%
% References:
%   [1] Di Renzo, M., Fu, L., and Urzay, J., HTR solver: An open-source
%       exascale-oriented task-based multi-GPU high-order code for
%       hypersonics aerothermodynamics, Comput. Phys. Commun, Vol. 255,
%       2020, p. 107262
%   [2] https://github.com/stanfordhpccenter/HTR-solver
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Apr 05 2023

clear; clc;

% Definitions
file_main = 'C:\Users\Alber\Mi unidad\Phd\Shocks\HTR\HIT\sample0'; % Path data parent folder
file_location_nodes = fullfile(file_main, 'cellCenter_grid\master.hdf'); % Path grid file
value = 0:5000:110000; % Folders tag
Tref = 300; % Temperature of reference [K]
mu_ref = 5.692750425533111 * 1e-4; % Dynamic viscosity of reference [kg/(m-s)] or [Pa-s]

% Miscellaneous
FLAG_SLICE = false; % Flag indicating to plot slice of the decomposed field
slice = 40;         % Slice to plot
FPS = 60;           % Frames per second

% Calculate decomposed velocity field and dissipation ratio
for i = length(value):-1:1
    fprintf('Case %2d:\n', i);
    file_location = fullfile(file_main, sprintf('fluid_iter%010d', value(i)), 'master.hdf');
    [eps_total(i), eps_ratio(i), time_norm(i),...
     K(i), K_ratio(i),...
     u_solenoidal(:, :, :, i), v_solenoidal(:, : , :, i),...
     u_compressive(:, :, :, i), v_compressive(:, : , :, i),...
     u_mean(i), v_mean(i), w_mean(i)] = helmholtz(file_location, file_location_nodes, Tref, mu_ref);
end

% Initialize figure
set_figure();
tiledlayout(2,2);

% Plot dissipation
ax1 = nexttile;
ax1 = set_figure(ax1);
plot(time_norm, eps_total, 'LineWidth', 1.2);
xlabel('$t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('$\epsilon$', 'Interpreter', 'latex');

% Plot dissipation ratio
ax2 = nexttile;
ax2 = set_figure(ax2);
plot(time_norm, eps_ratio, 'LineWidth', 1.2);
xlabel('$t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('$\epsilon_d / \epsilon_s$', 'Interpreter', 'latex');

% Plot turbulent kinetic energy
ax3 = nexttile;
ax3 = set_figure(ax3);
plot(time_norm, K, 'LineWidth', 1.2);
xlabel('$t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('$K$', 'Interpreter', 'latex');

% Plot turbulent kinetic energy ratio
ax4 = nexttile;
ax4 = set_figure(ax4);
plot(time_norm, K_ratio, 'LineWidth', 1.2);
xlabel('$t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('$K_d / K_s$', 'Interpreter', 'latex');

% Plot mean velocity field
ax5 = set_figure();
plot(time_norm, u_mean, time_norm, v_mean, time_norm, w_mean, 'LineWidth', 1.2);
xlabel('$t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('$\langle {u}_i \rangle$', 'Interpreter', 'latex');
leg_labels = {'$\langle u \rangle$', '$\langle v \rangle$', '$\langle w \rangle$'};
legend(ax5, leg_labels, 'Interpreter', 'latex');

% Plot slices of the decomposed field on the X-Y plane over time
if ~FLAG_SLICE
    return
end

set_figure();
for i = 1:length(value)
    plot_slice(u_solenoidal(:, : , :, i), v_solenoidal(:, : , :, i),...
               u_compressive(:, : , :, i), v_compressive(:, : , :, i),...
               size(u_compressive(:, :, :, 1)), slice);
    pause(1 / FPS);
end