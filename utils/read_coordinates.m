function [x, y, z, L] = read_coordinates(file_location_nodes)
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

    [coordinates_x, coordinates_y, coordinates_z,...
    sz_coordinates] = read_3D(file_location_nodes, 'centerCoordinates');

    x = reshape(coordinates_x(:, 1, 1), 1, sz_coordinates(1));
    y = reshape(coordinates_y(1, :, 1), 1, sz_coordinates(2));
    z = reshape(coordinates_z(1, 1, :), 1, sz_coordinates(3));

    % Get domain length
    L = max(x);
end