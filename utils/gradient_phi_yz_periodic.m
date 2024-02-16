function [u_dilatational, v_dilatational, w_dilatational] = gradient_phi_yz_periodic(phi, dx, dy, dz)
    % Compute the gradient assuming y and z directions are periodic
    
    % Initialization
    u_dilatational = zeros(size(phi));

    % Compute dudx (first and last points) - forward scheme
    u_dilatational(1, :, :) = (phi(2, :, :) - phi(1, :, :)) / dx;
    u_dilatational(end, :, :) = (phi(end, :, :) - phi(end-1, :, :)) / dx;
    % Compute dudx (inner points) - central scheme
    u_dilatational(2:end-1, :, :) = (phi(3:end, :, :) - phi(1:end-2, :, :)) / (2 * dx);
    
    % Compute dvdx - central scheme
    v_dilatational = (circshift(phi, [ 0, -1,  0]) - circshift(phi, [0, 1, 0])) ./ (2 * dy);

    % Compute dwdx - central scheme
    w_dilatational = (circshift(phi, [ 0,  0, -1]) - circshift(phi, [0, 0, 1])) ./ (2 * dz);
end