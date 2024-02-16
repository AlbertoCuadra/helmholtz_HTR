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