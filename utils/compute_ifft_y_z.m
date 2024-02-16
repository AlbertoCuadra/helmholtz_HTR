function ifft_result = compute_ifft_y_z(data, sz)

    % Initialize the result array
    ifft_result = zeros(sz);
    
    % Compute the IFFT along Y and Z dimensions for each X slice
    for i = 1:sz(1)
        ifft_result(i, :, :) = ifft2(data(i, :, :), "symmetric");  % 2D IFFT along Y and Z dimensions
    end

end