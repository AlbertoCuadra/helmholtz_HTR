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