function K = compute_tke_decomposed(u, v, w)
    % Compute Turbulent Kinetic Energy (TKE)
    %
    % Args:
    %     u (float): 3D array with the x-component of the velocity field
    %     v (float): 3D array with the y-component of the velocity field
    %     w (float): 3D array with the z-component of the velocity field
    %
    % Returns:
    %     K (float): 3D array with the TKE
    %
    % Example:
    %     K = compute_tke(u, v, w, rho);
    
    % K = 0.5 * mean(rho .* (u.^2 + v.^2 + w.^2), 'all');
    K = 0.5 * mean(u.^2 + v.^2 + w.^2, [2, 3]);
end