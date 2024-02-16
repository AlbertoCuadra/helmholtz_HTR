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