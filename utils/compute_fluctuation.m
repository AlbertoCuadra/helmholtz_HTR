function [fluctuation, average] = compute_fluctuation(value)
    % Compute fluctuation field and average from a given field
    %
    % Args:
    %     value (float): 3D array with a property field
    %
    % Returns:
    %     fluctuation (float): 3D array with the fluctuation of the property field
    %
    % Example:
    %     [fluctuation, average] = compute_fluctuation(u);

    average = mean(value, 'all');
    fluctuation = value - average;
end