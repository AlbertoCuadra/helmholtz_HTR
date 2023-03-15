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
% Last update Mar 15 2023

% Definitions
file_main = 'C:\Users\Alber\Mi unidad\Phd\Shocks\HTR\4605189\4605189\sample0'; % path data parent folder
file_location_nodes = fullfile(file_main, 'cellCenter_grid\master.hdf');   % path grid file
value = 52000:1000:63000; % folders tag

% Miscellaneous
slice = 40; % Slice to plot
FPS = 60;   % Frames per second

% Calculate decomposed velocity field and dissipation ratio 
for i = length(value):-1:1
    fprintf('Case %2d, ', i);
    file_location = fullfile(file_main, sprintf('fluid_iter%010d', value(i)), 'master.hdf'); % 0,0,0-79,79,79
    [eps_total(i), eps_ratio(i), time_norm(i), ...
     V_solenoidal_x(:, :, :, i), V_solenoidal_y(:, : , :, i),...
     V_compressive_x(:, :, :, i), V_compressive_y(:, : , :, i)...
    ] = helmholtz(file_location, file_location_nodes);
end

% Initialize figure
set_figure();
tiledlayout(1,2);

% Plot dissipation
ax1 = nexttile;
ax1 = set_figure(ax1);
plot(time_norm, eps_total);
xlabel('Time normalized, $t^* = t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('Dissipation, $\epsilon$', 'Interpreter', 'latex');
title(ax1, 'Dissipation', 'Interpreter', 'latex')

% Plot dissipation ratio
ax2 = nexttile;
ax2 = set_figure(ax2);
plot(time_norm, eps_ratio);
xlabel('Time normalized, $t^* = t / \epsilon_l$', 'Interpreter', 'latex');
ylabel('Dissipation ratio, $\epsilon_d / \epsilon_s$', 'Interpreter', 'latex');
title(ax2, 'Dissipation ratio', 'Interpreter', 'latex')

% Plot slices of the decomposed field on the X-Y plane over time
% set_figure();
% for i = 1:length(value)
%     plot_slice(V_solenoidal_x(:, : , :, i), V_solenoidal_y(:, : , :, i), V_compressive_x(:, : , :, i), V_compressive_y(:, : , :, i), size(V_compressive_x(:, :, :, 1)), slice);
%     pause(1 / FPS);
% end


