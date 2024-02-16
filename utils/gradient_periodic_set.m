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