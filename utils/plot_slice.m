function plot_slice(u_solenoidal, v_solenoidal, u_dilatational, v_dilatational, sz, slice)
    % Plot one slice of the decomposed field on X-Y plane
    %
    % Args:
    %     u_solenoidal (float): 3-D array with the solenoidal contribution of the velocity field in X-direction
    %     v_solenoidal (float): 3-D array with the solenoidal contribution of the velocity field in Y-direction
    %     u_dilatational (float): 3-D array with the dilatational contribution of the velocity field in X-direction
    %     v_dilatational (float): 3-D array with the dilatational contribution of the velocity field in Y-direction
    %     sz (3-element vector): The size of the input velocity field
    %     slice (int): The slice index to plot
    %
    % Example:
    %     plot_slice(u_solenoidal, v_solenoidal, u_dilatational, v_dilatational, sz, 1)

    % Set the quiver plot scaling factor
    plot_factor = 3.0;

    % Griddata scaling
    scaling_factor = 0.5;

    % Create a meshgrid for the X-Y plane
    [X, Y] = meshgrid(1:sz(1), 1:sz(2));
    
    % Transpose the matrices to swap Y and X dimensions
    X = permute(X, [2, 1]);
    Y = permute(Y, [2, 1]);

    % Griddata
    [XQ, YQ] = meshgrid(1:sz(1) * scaling_factor, 1:sz(2) * scaling_factor);
    XQ = permute(XQ, [2, 1]);
    YQ = permute(YQ, [2, 1]);
    
    u_dilatational_q = griddata(X, Y, u_dilatational(:, :, slice), XQ, YQ);
    v_dilatational_q = griddata(X, Y, v_dilatational(:, :, slice), XQ, YQ);
    u_solenoidal_q = griddata(X, Y, u_solenoidal(:, :, slice), XQ, YQ);
    v_solenoidal_q = griddata(X, Y, v_solenoidal(:, :, slice), XQ, YQ);
    
    % Set up the figure layout
    tiledlayout(2, 1);

    % Plot the solenoidal contribution
    ax1 = nexttile;
    ax1 = set_figure(ax1);
    quiver(XQ / scaling_factor, YQ / scaling_factor, u_solenoidal_q, v_solenoidal_q, 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax1, 'Solenoidal', 'Interpreter', 'latex')
    
    % Plot the dilatational contribution
    ax2 = nexttile;
    ax2 = set_figure(ax2);
    quiver(XQ / scaling_factor, YQ / scaling_factor, u_dilatational_q, v_dilatational_q, 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax2, 'Dilatational', 'Interpreter', 'latex')
end