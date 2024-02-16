function fft_result = compute_fft_y_z(data, sz)
    % Compute the FFT of the input data along Y and Z dimensions
    %
    % Args:
    %     data (float): 3D array with the input data
    %     sz (int): 3D array with the size of the input data
    %
    % Returns:
    %     fft_result (float): 3D array with the FFT of the input data along Y and Z dimensions
    %
    % Example:
    %     fft_result = compute_fft_y_z(data, sz)

    % Initialize the result array
    fft_result = zeros(sz);
    
    % Compute the FFT along Y and Z dimensions for each X slice
    for i = 1:sz(1)
        fft_result(i, :, :) = fft2(data(i, :, :));  % 2D FFT along Y and Z dimensions
    end

end