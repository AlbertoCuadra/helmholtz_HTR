function [eps_total, eps_ratio, time_norm, V_solenoidal_x,...
          V_solenoidal_y, V_compressive_x, V_compressive_y...
         ] = helmholtz(file_location, file_location_nodes)
    % Helmholtz-Hodge decomposition using fast Fourier transform [1].
    % 
    % This code is based on Ref. [2] and has been rewritten in MATLAB.
    %
    % Given the 3D velocity field, decompose into solenoidal and
    % compressive parts
    %
    % Note:
    %     1. Only for uniform grid with dx = dy = dz. Extensions are easy
    %        to implement.
    %     2. Helmholtz-Hodge decomposition performed with the spectral
    %        method should only apply to relatively smooth fields, i.e.,
    %        with little power on small scales.
    %     3. For even NX, NY, and NZ, decomposed fields can be complex,
    %        with the imaginary part coming from the real part of the kmode
    %        at Nyquist frequency. In principle, the Nyquist frequency
    %        kmode should be dropped when doing the first derivatives to
    %        maintain symmetry. See footnote on page 4 of [2]. However,
    %        when the field is smooth enough, the imaginary part caused by
    %        the Nyquist frequency kmode should be negligible.
    %
    % Args:
    %     file_location (string): Path to the data .hdf file
    %     file_location_nodes (string): Path to the grid .hdf file
    %
    % Returns:
    %     Tuple containing
    %
    %     * eps_total (float): Dissipation
    %     * eps_ratio (float): Dissipation ratio (solenoidal / compressive)
    %     * time_norm (float): Time normalized with the eddy turnover time
    %     * V_solenoidal_x (float): Solenoidal part of the velocity field in the x-axis
    %     * V_solenoidal_y (float): Solenoidal part of the velocity field in the y-axis
    %     * V_compressive_x (float): Compressive part of the velocity field in the x-axis
    %     * V_compressive_y (float): Compressive part of the velocity field in the y-axis
    %
    % References:
    %   [1] Johnson, S. G. (2011). Notes on FFT-based differentiation.
    %       MIT Applied Mathematics, Tech. Rep. 
    %       Available: http://math.mit.edu/~stevenj/fft-deriv.pdf
    %   [2] Xun Shi, Helmholtz-Hodge decomposition using fft (Python),
    %       Available: https://github.com/shixun22/helmholtz
    %
    %
    % @author: Alberto Cuadra Lara
    %          PhD Candidate - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                  
    % Last update Mar 15 2023
    
    % Definitions
    FPS = 60; % Frames per second
    mu = 1;   % Viscosity
    
    % Get coordinates
    [coordinates_x, coordinates_y, coordinates_z,...
     sz_coordinates] = read_3D(file_location_nodes, 'centerCoordinates');

    x = reshape(coordinates_x(:, 1, 1), 1, sz_coordinates(1));
    y = reshape(coordinates_y(1, :, 1), 1, sz_coordinates(2));
    z = reshape(coordinates_z(1, 1, :), 1, sz_coordinates(3));
    
    % Get domain length
    L = max(x);

    % Get velocity components
    [Vfx, Vfy, Vfz, sz] = read_3D(file_location, 'velocity');
    
    % Turbulent Kinetic Energy 
    TKE = 0.5 * mean(Vfx.^2 + Vfy.^2 + Vfz.^2, 'all');
    
    % Approximation of the integral length
    l = L / 5; 

    % Eddy turnover time
    eddy_turnover = l / (2 * TKE);

    % Time
    time = h5readatt(file_location, '/', 'simTime');

    % Time normalized with the Eddy turnover time
    time_norm = time / eddy_turnover;

    % Get density
    rho = read_data(file_location, 'rho');

    % Multiply by sqrt(rho)
    Vfx = sqrt(rho) .* Vfx;
    Vfy = sqrt(rho) .* Vfy;
    Vfz = sqrt(rho) .* Vfz;

    % Remove mean velocity
    Vfx = Vfx - mean(Vfx, 'all');
    Vfy = Vfy - mean(Vfy, 'all');
    Vfz = Vfz - mean(Vfz, 'all');
    
    % Get N-D fast Fourier transform (fft)
    vx_f = fftn(Vfx);
    vy_f = fftn(Vfy);
    vz_f = fftn(Vfz);

    % Get frequencies from the fft function
    kx = reshape(fftfreq(sz(1)), 1, 1, sz(1));
    ky = reshape(fftfreq(sz(2)), 1, sz(2), 1);
    kz = reshape(fftfreq(sz(3)), sz(3), 1, 1);

    % Get k^2
    k2 = kx.^2 + ky.^2 + kz.^2;

    % Avoid infinity value. We do not care about the k = 0 component
    k2(1, 1, 1) = 1;

    % Compute velocity divergence
    div_Vf_f = (vx_f .* kx + vy_f .* ky + vz_f .* kz);

    % Normalize
    V_compressive_overk = div_Vf_f ./ k2;

    % Get N-dimensional inverse discrete Fourier transform
    V_compressive_x = ifftn(V_compressive_overk .* kx);
    V_compressive_y = ifftn(V_compressive_overk .* ky);
    V_compressive_z = ifftn(V_compressive_overk .* kz);

    % Get solenoidal contributions
    V_solenoidal_x = Vfx - V_compressive_x;
    V_solenoidal_y = Vfy - V_compressive_y;
    V_solenoidal_z = Vfz - V_compressive_z;

    % Plot slices of the decomposed field on the X-Y plane
    % set_figure();
    % for slice = 1:sz(1)
    %     plot_slice(V_solenoidal_x, V_solenoidal_y, V_compressive_x, V_compressive_y, sz, slice);
    %     pause(1 / FPS);
    % end

    % Check if the solenoidal part is divergence-free
    divVs = ifftn((fftn(V_solenoidal_x) .* kx + ...
                   fftn(V_solenoidal_y) .* ky + ...
                   fftn(V_solenoidal_z) .* kz) * 1j * 2 * pi);
    fprintf('div_solenoidal max: %f\n', norm(divVs, 'fro'));
    
    % Compute velocity derivatives
    [dudx, dudy, dudz,...
     dvdx, dvdy, dvdz,...
     dwdx, dwdy, dwdz] = dxdy_periodic_3D_set(Vfx, Vfy, Vfz, x, y, z);
    
    % Compute divergence velocity | div(V(x, y, z)) = dudx + dvdy + dwdz
    div_v = dudx + dvdy + dwdz;

    % Calculate shear-stress tensor tau_ij obeying Newton's relation and
    % Stokes' hypothesis
    tau_11 = mu * (2 * dudx);
    tau_22 = mu * (2 * dvdy);
    tau_33 = mu * (2 * dwdz);
    tau_12 = mu * (dudy + dvdx - 2/3 * div_v);
    tau_13 = mu * (dudz + dwdx - 2/3 * div_v);
    tau_23 = mu * (dvdz + dwdy - 2/3 * div_v);

    % Calculate dissipation components eps_alpha = avg(tau_ij * w_i,alpha / sqrt(rho))
    eps_solenoidal = mean(dissipation(V_solenoidal_x, V_solenoidal_y, V_solenoidal_z), 'all');
    eps_compressive = mean(dissipation(V_compressive_x, V_compressive_y, V_compressive_z), 'all');

    % Compute dissipation ratio eps_compressive / eps_solenoidal
    eps_ratio = eps_compressive ./ eps_solenoidal;
    
    % Compute dissipation
    eps_total =  mean(dissipation(Vfx, Vfy, Vfz), 'all');

    % NESTED FUNCTIONS
    function value = dissipation(Vfx, Vfy, Vfz)
        [dudx, dudy, dudz,...
         dvdx, dvdy, dvdz,...
         dwdx, dwdy, dwdz] = dxdy_periodic_3D_set(Vfx, Vfy, Vfz, x, y, z);

        value = ...
            (tau_11 .* dudx + tau_12 .* dudy + tau_13 .* dudz + ...
             tau_12 .* dvdx + tau_22 .* dvdy + tau_23 .* dvdz + ...
             tau_13 .* dwdx + tau_23 .* dwdy + tau_33 .* dwdz) ./ sqrt(rho);
    end

end

% SUB-PASS FUNCTIONS
function value = read_data(file_location, properties)
    % Read data
    value = h5read(file_location, ['/', properties]);
    % Reshape
    sz = size(value);
    N_min = min(sz);
    if sz(1) ~= N_min
        value = value';
    end

end

function [x, y, z, sz] = read_3D(file_location, properties)
    % Read velocity data
    value = h5read(file_location, ['/', properties]);
    % Reshape
    sz = size(value);
    value = reshape(value, sz);
    % Get Dimensions
    sz = sz(2:end);
    % Get components
    x(:, :, :) = value(1, :, :, :);
    y(:, :, :) = value(2, :, :, :);
    z(:, :, :) = value(3, :, :, :);
end

function dxdy = dxdy_periodic(x, y)
    % Compute first central derivate using an uniform grid and considering
    % a periodic box
    %
    % Args:
    %     x (float): Grid values
    %     y (float): Values for the corresponding grid
    %
    % Returns:
    %     dxdy (float): Value of the first derivate for the given grid and its corresponding values
    
    % Definitions
    n = length(y); 
    % Get step
    h = y(2) - y(1);
    % Preallocate derivative vector
    dxdy = zeros(1, n);
    % Calculate first and last points considering a periodic box
    dxdy(1) = (x(2) - x(end)) / (2 * h);
    dxdy(end) = (x(1) - x(end - 1)) / (2 * h);
    % Calculate the remaining points
    for i = 2:n - 1
        dxdy(i) = (x(i + 1) - x(i - 1)) / (2 * h);
    end

end

function dxdy_3D = dxdy_periodic_3D(x, y, direction)
    % Compute first central derivate using an uniform grid of a 3D field
    % considering a periodic box
    %
    % Args:
    %     x (float): Grid values
    %     y (float): Values for the corresponding grid
    %
    % Returns:
    %     dxdy (float): Value of the first derivate for the given grid and its corresponding values
    
    % Definitions
    sz = size(x);
    % Calculations
    switch direction
        case 1
            for k = sz(3):-1:1
                for j = sz(2):-1:1
                    dxdy_3D(:, j, k) = dxdy_periodic(x(:, j, k), y);
                end
            end

        case 2
            for k = sz(3):-1:1
                for i = sz(1):-1:1
                    dxdy_3D(i, :, k) = dxdy_periodic(x(i, :, k), y);
                end
            end

        case 3
            for i = sz(1):-1:1
                for j = sz(2):-1:1
                    dxdy_3D(i, j, :) = dxdy_periodic(x(i, j, :), y);
                end
            end
    end
end

function [dudx, dudy, dudz,...
          dvdx, dvdy, dvdz,...
          dwdx, dwdy, dwdz] = dxdy_periodic_3D_set(Vfx, Vfy, Vfz, x, y, z)
    % Compute first derivative of a 3D field over the three axis
    % considering a periodic box

    % 1. dudx_i
    dudx = dxdy_periodic_3D(Vfx, x, 1);
    dudy = dxdy_periodic_3D(Vfx, y, 2);
    dudz = dxdy_periodic_3D(Vfx, z, 3);

    % 2. dvdx_i
    dvdx = dxdy_periodic_3D(Vfy, x, 1);
    dvdy = dxdy_periodic_3D(Vfy, y, 2);
    dvdz = dxdy_periodic_3D(Vfy, z, 3);

    % 3. dwdx_i
    dwdx = dxdy_periodic_3D(Vfz, x, 1);
    dwdy = dxdy_periodic_3D(Vfz, y, 2);
    dwdz = dxdy_periodic_3D(Vfz, z, 3);
end