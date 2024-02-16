function [eps_total, eps_ratio, time_norm,...
          K, K_ratio, u_solenoidal, v_solenoidal,...
          u_dilatational, v_dilatational,...
          u_mean, v_mean, w_mean,...
          K_dilatational, K_solenoidal, chi] = helmholtz(file_location,...
                                                file_location_nodes,...
                                                T_ref,...
                                                mu_ref)

    % Helmholtz-Hodge decomposition of a 3D velocity field into its
    % solenoidal and dilatational parts using fast Fourier transform [1].
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
    %     * eps_ratio (float): Dissipation ratio (dilatational / solenoidal)
    %     * time_norm (float): Time normalized with the eddy turnover time
    %     * u_solenoidal (float): Solenoidal part of the velocity field in the x-axis
    %     * v_solenoidal (float): Solenoidal part of the velocity field in the y-axis
    %     * u_dilatational (float): Dilatational part of the velocity field in the x-axis
    %     * v_dilatational (float): Dilatational part of the velocity field in the y-axis
    %     * u_mean (float): Mean of the velocity field in the x-axis
    %     * v_mean (float): Mean of the velocity field in the y-axis
    %     * w_mean (float): Mean of the velocity field in the z-axis
    %     * K_dilatational (float): Dilatational contribution of the TKE
    %     * K_solenoidal (float): Solenoidal contribution of the TKE
    %     * chi (float): Ratio density deviations / velocity perturbations
    %
    % Example:
    %     [eps_total, eps_ratio, time_norm,...
    %      K, K_ratio, u_solenoidal, v_solenoidal,...
    %      u_dilatational, v_dilatational,...
    %      u_mean, v_mean, w_mean,...
    %      K_dilatational, K_solenoidal, chi] = helmholtz(file_location,...
    %                                            file_location_nodes,...
    %                                            T_ref,...
    %                                            mu_ref)
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
    %          Postdoctoral researcher - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                  
    % Last update Feb 16 2024

    % Definitions
    FPS = 60;  % Frames per second
    S = 110.4; % Sutherland constant [K]
    gamma_mean = 1.4; % Adiabatic index

    % Get coordinates
    [x, y, z, L] = read_coordinates(file_location_nodes);

    % Get density, temperature, and pressure
    rho = read_data(file_location, 'rho');
    T = read_data(file_location, 'temperature');
    p = read_data(file_location, 'pressure');
    [rho_delta, rho_mean] = compute_fluctuation(rho);
    [T_delta, T_mean] = compute_fluctuation(T);
    [p_delta, p_mean] = compute_fluctuation(p);

    % Compute mean sound speed
    sound_mean = sqrt(gamma_mean * p_mean ./ rho_mean);

    % Compute dynamic viscosity field using the Sutherland's law
    mu = compute_mu_sutherland(T, S, T_ref, mu_ref);

    % Get velocity components
    [wx, wy, wz, sz] = read_3D(file_location, 'velocity');

    % Get velocity components * sqrt(rho)
    u = sqrt(rho) .* wx;
    v = sqrt(rho) .* wy;
    w = sqrt(rho) .* wz;

    % Get velocity fluctuations
    [wx, wx_mean] = compute_fluctuation(wx);
    [wy, wy_mean] = compute_fluctuation(wy);
    [wz, wz_mean] = compute_fluctuation(wz);
    [u, u_mean] = compute_fluctuation(u);
    [v, v_mean] = compute_fluctuation(v);
    [w, w_mean] = compute_fluctuation(w);

    % Module and root mean square of the velocity field * sqrt(rho)
    vel = sqrt(wx.^2 + wy.^2 + wz.^2);
    vel_rms = sqrt( (wx.^2 + wy.^2 + wz.^2) / 3);

    % Get turbulent Mach number
    Mt = sqrt(3) * vel_rms / sound_mean;
    Mt_mean = mean(Mt, 'all');

    fprintf('Turbulent Mach:        %.2e\n', Mt_mean);

    % Compute ratio density deviations / velocity perturbations (for LIA)
    % chi_ast = mean(rho_delta .* vel, 'all') / (rho_mean * sound_mean);
    % chi = chi_ast / Mt_mean^2;
    % chi = -abs(mean(rho_delta ./ Mt, 'all') / rho_mean);

    % vel = u;
    chi = (mean(rho_delta .* vel, 'all') / (rho_mean * sound_mean)) / (mean(vel .* vel, 'all') / sound_mean^2);
    % chi = -(mean(T_delta .* vel, 'all') / (T_mean * sound_mean)) / (mean(vel .* vel, 'all') / sound_mean^2);

    fprintf('chi (Cuadra):          %.2e\n', chi);
    fprintf('chi (Sinha):           %.2e\n', chi * (wx_mean / sound_mean));

    % Approximation of the integral length
    l = L / 6; 

    % Eddy turnover time
    eddy_turnover = l / mean(vel_rms, 'all');

    % Time
    time = h5readatt(file_location, '/', 'simTime');

    % Time normalized with the Eddy turnover time
    time_norm = time / eddy_turnover;

    % Get N-D fast Fourier transform (fft)
    U = fftn(u);
    V = fftn(v);
    W = fftn(w);

    % Get wave numbers
    kx = fftfreq(sz(1));
    ky = fftfreq(sz(2));
    kz = fftfreq(sz(3));

    % Get grid wave numbers
    [KX, KY, KZ] = ndgrid(kx, ky, kz);

    % Get k^2
    K2 = KX.^2 + KY.^2 + KZ.^2;

    % Avoid infinity value. We do not care about the k = 0 component
    K2(1, 1, 1) = 1;

    % Compute velocity divergence
    div = (U .* KX + V .* KY + W .* KZ);

    % Compute the Helmholtz decomposition (dilatational)
    H = div ./ K2;

    % Get dilatational contributions (curl-free)
    u_dilatational = ifftn(H .* KX);
    v_dilatational = ifftn(H .* KY);
    w_dilatational = ifftn(H .* KZ);

    % Get solenoidal contributions (divergence-free)
    u_solenoidal = u - u_dilatational;
    v_solenoidal = v - v_dilatational;
    w_solenoidal = w - w_dilatational;

    % Plot slices of the decomposed field on the X-Y plane
    % set_figure();
    % for slice = 1:sz(1)
    %     plot_slice(u_solenoidal, v_solenoidal, u_dilatational, v_dilatational, sz, slice);
    %     pause(1 / FPS);
    % end

    % Check if the solenoidal part is divergence-free
    div_solenoidal = ifftn((fftn(u_solenoidal) .* KX + ...
                            fftn(v_solenoidal) .* KY + ...
                            fftn(w_solenoidal) .* KZ) * 1i * 2 * pi);

    fprintf('div_solenoidal max:    %.2e\n', max(abs(div_solenoidal(:))));

    % Check if the dilatational part is curl-free
    curl_dilatational = ifftn((fftn(w_dilatational) .* KY - fftn(v_dilatational) .* KZ + ...
                               fftn(u_dilatational) .* KZ - fftn(w_dilatational) .* KX + ...
                               fftn(v_dilatational) .* KX - fftn(u_dilatational) .* KY) * 1i * 2 * pi);

    fprintf('curl_dilatational max: %.2e\n', max(abs(curl_dilatational(:))));

    % Compute velocity derivatives
    [dudx, dudy, dudz,...
    dvdx, dvdy, dvdz,...
    dwdx, dwdy, dwdz] = gradient_periodic_set(u./sqrt(rho), v./sqrt(rho), w./sqrt(rho), x, y, z);

    % Compute divergence velocity | div( vf(x, y, z) ) = dudx + dvdy + dwdz
    div_vf = dudx + dvdy + dwdz;

    % Calculate shear-stress tensor tau_ij assuming Newton's relation
    % and Stokes' hypothesis
    tau.tau_11 = mu .* (2 * dudx);
    tau.tau_22 = mu .* (2 * dvdy);
    tau.tau_33 = mu .* (2 * dwdz);
    tau.tau_12 = mu .* (dudy + dvdx - 2/3 * div_vf);
    tau.tau_13 = mu .* (dudz + dwdx - 2/3 * div_vf);
    tau.tau_23 = mu .* (dvdz + dwdy - 2/3 * div_vf);

    % Calculate pressure dilatation
    p_dilatational = mean(p_delta .* div_vf, 'all');

    % Calculate dissipation components eps_alpha = avg(tau_ij * w_i,alpha / sqrt(rho))
    eps_solenoidal = mean(dissipation(u_solenoidal, v_solenoidal, w_solenoidal, x, y, z, tau), 'all');
    eps_dilatational = mean(dissipation(u_dilatational, v_dilatational, w_dilatational, x, y, z, tau), 'all');

    % Compute dissipation ratio eps_dilatational / eps_solenoidal
    eps_ratio = eps_dilatational ./ eps_solenoidal;

    % Compute dissipation
    eps_total =  mean(dissipation(u, v, w, x, y, z, tau), 'all');

    fprintf('Total dissipation:     %.2e\n', eps_total);

    % Compute dilatational and solenoidal contributions of the TKE
    K_solenoidal = compute_tke(u_solenoidal, v_solenoidal, w_solenoidal, rho);
    K_dilatational = compute_tke(u_dilatational, v_dilatational, w_dilatational, rho);

    % Turbulent Kinetic Energy 
    K = compute_tke(u, v, w, rho);

    % Ratio dilatational to total TKE
    ratio_dilatational = K_dilatational / K;

    fprintf('Ratio dilatational:    %.2e\n\n', ratio_dilatational);

    % Compute TKE ratio
    K_ratio = K_dilatational ./ K_solenoidal;

    % Check if system is in equilibrium
    eps_dilatational_equil = p_dilatational + 2 * sound_mean * K_dilatational;
end