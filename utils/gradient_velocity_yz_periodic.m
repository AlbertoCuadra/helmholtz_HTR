function [grad_u, grad_v, grad_w] = gradient_velocity_yz_periodic(u, v, w, dx, dy, dz)
    % Compute the gradient assuming y and z directions are periodic
    
    % Initialization
    grad_u = zeros(size(u));

    % Compute dudx (first and last points) - forward scheme
    grad_u(1, :, :) = (u(2, :, :) - u(1, :, :)) / dx;
    grad_u(end, :, :) = (u(end, :, :) - u(end-1, :, :)) / dx;

    % Compute dudx (inner points) - central scheme
    grad_u(2:end-1, :, :) = (u(3:end, :, :) - u(1:end-2, :, :)) / (2 * dx);
    
    % Compute dvdx - central scheme
    grad_v = (circshift(v, [ 0, -1,  0]) - circshift(v, [0, 1, 0])) ./ (2 * dy);

    % Compute dwdx - central scheme
    grad_w = (circshift(w, [ 0,  0, -1]) - circshift(w, [0, 0, 1])) ./ (2 * dz);
end