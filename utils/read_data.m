function value = read_data(file_location, property)
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

    value = h5read(file_location, ['/', property]);
    % Reshape
    sz = size(value);
    N_min = min(sz);
    if sz(1) ~= N_min
        value = value';
    end

end