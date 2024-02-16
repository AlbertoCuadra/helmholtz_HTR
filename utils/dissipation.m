function value = dissipation(u, v, w, x, y, z, tau)
    % Compute dissipation
    %
    % Args:
    %   u (float): 3D array with the x-component of the velocity field
    %   v (float): 3D array with the y-component of the velocity field
    %   w (float): 3D array with the z-component of the velocity field
    %   x (float): 3D array with the x-coordinates
    %   y (float): 3D array with the y-coordinates
    %   z (float): 3D array with the z-coordinates
    %   tau (struct): Structure with the 9 components of the shear-stress tensor
    %
    % Returns:
    %   value (float): 3D array with the dissipation
    %
    % Example:
    %   value = dissipation(u, v, w, x, y, z, tau);
    
    % Compute velocity derivatives
    [dudx, dudy, dudz,...
     dvdx, dvdy, dvdz,...
     dwdx, dwdy, dwdz] = gradient_periodic_set(u, v, w, x, y, z);
    
    % Compute dissipation
    value = ...
        ( (tau.tau_11 .* dudx + tau.tau_12 .* dudy + tau.tau_13 .* dudz) + ...
          (tau.tau_12 .* dvdx + tau.tau_22 .* dvdy + tau.tau_23 .* dvdz) + ...
          (tau.tau_13 .* dwdx + tau.tau_23 .* dwdy + tau.tau_33 .* dwdz) );
end