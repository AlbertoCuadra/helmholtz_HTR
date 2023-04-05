function [eps_total, eps_ratio, time_norm,...
          K, K_ratio, u_solenoidal, v_solenoidal,...
          u_compressive, v_compressive,...
          u_mean, v_mean, w_mean] = helmholtz(file_location,...
                                              file_location_nodes,...
                                              T_ref,...
                                              mu_ref)
    % Helmholtz-Hodge decomposition of a 3D velocity field into its
    % solenoidal and compressive parts using fast Fourier transform [1].
    % 
    % The decomposition is performed with the spectral method, which is
    % only suitable for relatively smooth fields, i.e., with little power
    % on small scales. The code assumes that the grid is uniform with
    % dx = dy = dz.
    %
    % This code is based on Ref. [2] and has been rewritten in MATLAB.
    %
    % Notes:
    %     For even NX, NY, and NZ, decomposed fields can be complex,
    %     with the imaginary part coming from the real part of the kmode
    %     at Nyquist frequency. In principle, the Nyquist frequency
    %     kmode should be dropped when doing the first derivatives to
    %     maintain symmetry. See footnote on page 4 of [2]. However,
    %     when the field is smooth enough, the imaginary part caused by
    %     the Nyquist frequency kmode should be negligible.
    %
    % Args:
    %     file_location (char): Path to the data .hdf file
    %     file_location_nodes (char): Path to the grid .hdf file
    %     T_ref (float): Temperature of reference [K]
    %     mu_ref (float): Dynamic viscosity of reference [kg/(m-s)] or [Pa-s]
    %
    % Returns:
    %     Tuple containing
    %
    %     * eps_total (float): Total dissipation
    %     * eps_ratio (float): Dissipation ratio (solenoidal / compressive)
    %     * time_norm (float): Time normalized with the eddy turnover time
    %     * u_solenoidal (float): Solenoidal part of the velocity field in the x-axis
    %     * v_solenoidal (float): Solenoidal part of the velocity field in the y-axis
    %     * u_compressive (float): Compressive part of the velocity field in the x-axis
    %     * v_compressive (float): Compressive part of the velocity field in the y-axis
    %     * u_mean (float): Mean of the velocity field in the x-axis
    %     * v_mean (float): Mean of the velocity field in the y-axis
    %     * w_mean (float): Mean of the velocity field in the z-axis
    %
    % Examples:
    %     [eps_total, eps_ratio, time_norm,...
    %      K, K_ratio, u_solenoidal, v_solenoidal,...
    %      u_compressive, v_compressive,...
    %      u_mean, v_mean, w_mean] = helmholtz(file_location,...
    %                                          file_location_nodes,...
    %                                          T_ref,...
    %                                          mu_ref)
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
    % Last update Apr 05 2023
    
    % Definitions
    FPS = 60; % Frames per second
    S = 110.4; % Sutherland constant [K]

    % Get coordinates
    [coordinates_x, coordinates_y, coordinates_z,...
     sz_coordinates] = read_3D(file_location_nodes, 'centerCoordinates');

    x = reshape(coordinates_x(:, 1, 1), 1, sz_coordinates(1));
    y = reshape(coordinates_y(1, :, 1), 1, sz_coordinates(2));
    z = reshape(coordinates_z(1, 1, :), 1, sz_coordinates(3));
    
    % Get domain length
    L = max(x);

    % Get velocity components
    [u, v, w, sz] = read_3D(file_location, 'velocity');
    
    % Get mean velocity field
    u_mean = mean(u, "all");
    v_mean = mean(v, "all");
    w_mean = mean(w, "all");
    
    % Remove mean velocity
    u = u - mean(u, 'all');
    v = v - mean(v, 'all');
    w = w - mean(w, 'all');

    % Get density
    rho = read_data(file_location, 'rho');

    % Turbulent Kinetic Energy 
    K = compute_tke(u, v, w, rho);
    
    % Approximation of the integral length
    l = L / 6; 
    
    % Velocity rms
    vel_rms = sqrt(mean(u.^2 + v.^2 + w.^2, "all") / 3);

    % Eddy turnover time
    eddy_turnover = l / vel_rms;

    % Time
    time = h5readatt(file_location, '/', 'simTime');

    % Time normalized with the Eddy turnover time
    time_norm = time / eddy_turnover;

    % Multiply by sqrt(rho)
    u = sqrt(rho) .* u;
    v = sqrt(rho) .* v;
    w = sqrt(rho) .* w;
    
    % Get N-D fast Fourier transform (fft)
    U = fftn(u);
    V = fftn(v);
    W = fftn(w);

    % Get wave numbers
    kx = fftfreq(sz(1));
    ky = fftfreq(sz(2));
    kz = fftfreq(sz(3));

    % Get grid wave numbers
    [KX, KY, KZ] = meshgrid(kx, ky, kz);
    
    % Get k^2
    K2 = KX.^2 + KY.^2 + KZ.^2;

    % Avoid infinity value. We do not care about the k = 0 component
    K2(1, 1, 1) = 1;

    % Compute velocity divergence
    div = (U .* KX + V .* KY + W .* KZ);
    
    % Compute the Helmholtz decomposition (compressive)
    H = div ./ K2;

    % Get compressive contributions (curl-free)
    u_compressive = ifftn(H .* KX);
    v_compressive = ifftn(H .* KY);
    w_compressive = ifftn(H .* KZ);

    % Get solenoidal contributions (divergence-free)
    u_solenoidal = u - u_compressive;
    v_solenoidal = v - v_compressive;
    w_solenoidal = w - w_compressive;

    % Plot slices of the decomposed field on the X-Y plane
    % set_figure();
    % for slice = 1:sz(1)
    %     plot_slice(u_solenoidal, v_solenoidal, u_compressive, v_compressive, sz, slice);
    %     pause(1 / FPS);
    % end

    % Check if the solenoidal part is divergence-free
    div_solenoidal = ifftn((fftn(u_solenoidal) .* KX + ...
                            fftn(v_solenoidal) .* KY + ...
                            fftn(w_solenoidal) .* KZ) * 1i * 2 * pi);
    
    fprintf('div_solenoidal max:   %.2e\n', max(abs(div_solenoidal(:))));

    % Check if the compressive part is curl-free
    curl_compressive = ifftn((fftn(w_compressive) .* KY - fftn(v_compressive) .* KZ + ...
                              fftn(u_compressive) .* KZ - fftn(w_compressive) .* KX + ...
                              fftn(v_compressive) .* KX - fftn(u_compressive) .* KY) * 1i * 2 * pi);
    
    fprintf('curl_compressive max: %.2e\n', max(abs(curl_compressive(:))));

    % Compute velocity derivatives
    [dudx, dudy, dudz,...
     dvdx, dvdy, dvdz,...
     dwdx, dwdy, dwdz] = gradient_periodic_set(u, v, w, x, y, z);
    
    % Compute divergence velocity | div( vf(x, y, z) ) = dudx + dvdy + dwdz
    div_vf = dudx + dvdy + dwdz;
    
    % Get temperature field
    T = read_data(file_location, 'temperature');
    
    % Compute dynamic viscosity field using the Sutherland's law
    mu = compute_mu_sutherland(T, S, T_ref, mu_ref);

    % Calculate shear-stress tensor tau_ij assuming Newton's relation
    % and Stokes' hypothesis
    tau_11 = mu .* (2 * dudx);
    tau_22 = mu .* (2 * dvdy);
    tau_33 = mu .* (2 * dwdz);
    tau_12 = mu .* (dudy + dvdx - 2/3 * div_vf);
    tau_13 = mu .* (dudz + dwdx - 2/3 * div_vf);
    tau_23 = mu .* (dvdz + dwdy - 2/3 * div_vf);

    % Calculate dissipation components eps_alpha = avg(tau_ij * w_i,alpha / sqrt(rho))
    eps_solenoidal = mean(dissipation(u_solenoidal, v_solenoidal, w_solenoidal), 'all');
    eps_compressive = mean(dissipation(u_compressive, v_compressive, w_compressive), 'all');

    % Compute dissipation ratio eps_compressive / eps_solenoidal
    eps_ratio = eps_compressive ./ eps_solenoidal;
    
    % Compute dissipation
    eps_total =  mean(dissipation(u, v, w), 'all');

    fprintf('Total dissipation:    %.2e\n\n', eps_total);
    
    % Compute compressive and solenoidal contributions of the TKE
    K_solenoidal = compute_tke(u_solenoidal, v_solenoidal, w_solenoidal, 1);
    K_compressive = compute_tke(u_compressive, v_compressive, w_compressive, 1);
    
    % Compute TKE ratio
    K_ratio = K_solenoidal ./ K_compressive;
    
    % NESTED FUNCTIONS
    function value = dissipation(u, v, w)
        % Compute dissipation
        %
        % Args:
        %   u (float): 3D array with the x-component of the velocity field
        %   v (float): 3D array with the y-component of the velocity field
        %   w (float): 3D array with the z-component of the velocity field
        %
        % Returns:
        %   value (float):  3D array with the dissipation
        %
        % Example:
        %   value = dissipation(u, v, w);

        [dudx, dudy, dudz,...
         dvdx, dvdy, dvdz,...
         dwdx, dwdy, dwdz] = gradient_periodic_set(u, v, w, x, y, z);

        value = ...
            (tau_11 .* dudx + tau_12 .* dudy + tau_13 .* dudz + ...
             tau_12 .* dvdx + tau_22 .* dvdy + tau_23 .* dvdz + ...
             tau_13 .* dwdx + tau_23 .* dwdy + tau_33 .* dwdz) ./ sqrt(rho);
    end

end

% SUB-PASS FUNCTIONS
function value = read_data(file_location, property)
    % Read data from a .hdf file
    %
    % Args:
    %     file_location (char): Path to the .hdf file
    %     property (char): Name of the property to read
    %
    % Returns:
    %     value (float): 3D array with the property field
    %
    % Example:
    %     value = read_data(file_location, 'rho');

    value = h5read(file_location, ['/', property]);
    % Reshape
    sz = size(value);
    N_min = min(sz);
    if sz(1) ~= N_min
        value = value';
    end

end

function [x, y, z, sz] = read_3D(file_location, property)
    % Read velocity data
    %
    % Args:
    %     file_location (char): Path to the .hdf file
    %     property (char): Name of the property to read
    %
    % Returns:
    %     x (float): 3D array with the x-component of the property field
    %     y (float): 3D array with the y-component of the property field
    %     z (float): 3D array with the z-component of the property field
    %     sz (float): Size of the 3D array
    %
    % Example:
    %     [u, v, w, sz] = read_3D(file_location, 'velocity');
    
    value = h5read(file_location, ['/', property]);
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

function mu = compute_mu_sutherland(T, S, T_ref, mu_ref)
    % Compute dynamic viscosity using the Sutherland's law
    %
    % Args:
    %     T (float): 3D array with the temperature field
    %     S (float): Sutherland's constant [K]
    %     T_ref (float): Temperature of reference [K]
    %     mu_ref (float): Dynamic viscosity of reference [kg/(m-s)] or [Pa-s]
    %
    % Returns:
    %     mu (float): 3D array with the dynamic viscosity
    %
    % Example:
    %     mu = compute_mu_sutherland(T, S, T_ref, mu_ref);

    C1 = mu_ref * ((273.15 + 110.4) / T_ref / (273.15 / T_ref)^(3/2));
    mu =  C1 * T.^(3/2) ./ (T + S / T_ref);
end

function K = compute_tke(u, v, w, rho)
    % Compute Turbulent Kinetic Energy (TKE)
    %
    % Args:
    %     u (float): 3D array with the x-component of the velocity field
    %     v (float): 3D array with the y-component of the velocity field
    %     w (float): 3D array with the z-component of the velocity field
    %     rho (float): 3D array with the density field
    %
    % Returns:
    %     K (float): 3D array with the TKE
    %
    % Example:
    %     K = compute_tke(u, v, w, rho);
    
    K = 0.5 * mean(rho .* (u.^2 + v.^2 + w.^2), 'all');
end

function [dudx, dudy, dudz,...
          dvdx, dvdy, dvdz,...
          dwdx, dwdy, dwdz] = gradient_periodic_set(u, v, w, x, y, z)
    % Compute first derivative of a 3D field over the three axis
    % considering a periodic box
    %
    % Args:
    %     u (float): 3D array with the x-component of the velocity field
    %     v (float): 3D array with the y-component of the velocity field
    %     w (float): 3D array with the z-component of the velocity field
    %     x (float): 3D array with the x-coordinate of the grid
    %     y (float): 3D array with the y-coordinate of the grid
    %     z (float): 3D array with the z-coordinate of the grid
    %
    % Returns:
    %     dudx (float): 3D array with the first derivative of u over x
    %     dudy (float): 3D array with the first derivative of u over y
    %     dudz (float): 3D array with the first derivative of u over z
    %     dvdx (float): 3D array with the first derivative of v over x
    %     dvdy (float): 3D array with the first derivative of v over y
    %     dvdz (float): 3D array with the first derivative of v over z
    %     dwdx (float): 3D array with the first derivative of w over x
    %     dwdy (float): 3D array with the first derivative of w over y
    %     dwdz (float): 3D array with the first derivative of w over z
    %
    % Example:
    %     [dudx, dudy, dudz,...
    %      dvdx, dvdy, dvdz,...
    %      dwdx, dwdy, dwdz] = gradient_periodic_set(u, v, w, x, y, z);
    
    % 1. dudx_i
    [dudx, dudy, dudz] = gradient_periodic(u, x, y, z);
    % 2. dvdx_i
    [dvdx, dvdy, dvdz] = gradient_periodic(v, x, y, z);
    % 3. dwdx_i
    [dwdx, dwdy, dwdz] = gradient_periodic(w, x, y, z);
end

function [grad_x, grad_y, grad_z] = gradient_periodic(f, x, y, z)
    % Computes the gradient of a three-dimensional scalar field using
    % second-order central finite differences considering a periodic box
    %
    % Args:
    %     f (float): Scalar field to compute the gradient of
    %     x (float): 3D array with the x-coordinate of the grid
    %     y (float): 3D array with the y-coordinate of the grid
    %     z (float): 3D array with the z-coordinate of the grid
    %
    % Returns:
    %     grad_x (float): 3D array with the gradient of f over x
    %     grad_y (float): 3D array with the gradient of f over y
    %     grad_z (float): 3D array with the gradient of f over z
    %
    % Example:
    %     [grad_x, grad_y, grad_z] = gradient_periodic(f, x, y, z);
    
    % Get grid spacing in the three directions
    hx = x(2) - x(1);
    hy = y(2) - y(1);
    hz = z(2) - z(1);

    % Compute the gradient in each direction using central finite differences
    % considering periodic boundary conditions
    grad_x = (circshift(f, [-1,  0,  0]) - circshift(f, [1, 0, 0])) ./ (2 * hx);
    grad_y = (circshift(f, [ 0, -1,  0]) - circshift(f, [0, 1, 0])) ./ (2 * hy);
    grad_z = (circshift(f, [ 0,  0, -1]) - circshift(f, [0, 0, 1])) ./ (2 * hz);
end