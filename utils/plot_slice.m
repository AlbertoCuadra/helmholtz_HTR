function plot_slice(u_solenoidal, v_solenoidal, u_compressive, v_compressive, sz, slice)
    % Plot one slice of the decomposed field on X-Y plane
    %
    % Args:
    %     u_solenoidal (3-D array): The solenoidal contribution of the velocity field in X-direction
    %     v_solenoidal (3-D array): The solenoidal contribution of the velocity field in Y-direction
    %     u_compressive (3-D array): The compressive contribution of the velocity field in X-direction
    %     v_compressive (3-D array): The compressive contribution of the velocity field in Y-direction
    %     sz (3-element vector): The size of the input velocity field
    %     slice (int): The slice index to plot
    %
    % Returns:
    %     None

    % Set the quiver plot scaling factor
    plot_factor = 3;

    % Create a meshgrid for the X-Y plane
    [X, Y] = meshgrid(1:sz(1), 1:sz(2));

    % Set up the figure layout
    tiledlayout(1,2);

    % Plot the solenoidal contribution
    ax1 = nexttile;
    ax1 = set_figure(ax1);
    quiver(X, Y, u_solenoidal(:, :, slice), v_solenoidal(:, :, slice), 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax1, 'Solenoidal', 'Interpreter', 'latex')
    
    % Plot the compressive contribution
    ax2 = nexttile;
    ax2 = set_figure(ax2);
    quiver(X, Y, u_compressive(:, :, slice), v_compressive(:, :, slice), 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax2, 'Compressive', 'Interpreter', 'latex')
end