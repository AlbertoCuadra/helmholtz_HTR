function plot_slice(V_solenoidal_x, V_solenoidal_y, V_compressive_x, V_compressive_y, sz, slice)
    % Plot one slice of the decomposed field on X-Y plane
    plot_factor = 3;
    [X, Y] = meshgrid(1:sz(1), 1:sz(2));

    tiledlayout(1,2);

    ax1 = nexttile;
    ax1 = set_figure(ax1);
    quiver(X, Y, V_solenoidal_x(:, :, slice), V_solenoidal_y(:, :, slice), 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax1, 'Solenoidal', 'Interpreter', 'latex')
    
    ax2 = nexttile;
    ax2 = set_figure(ax2);
    quiver(X, Y, V_compressive_x(:, :, slice), V_compressive_y(:, :, slice), 'AutoScale', 'on', 'AutoScaleFactor', plot_factor);
    xlim([1, sz(1)]); ylim([1, sz(2)]);
    title(ax2, 'Compressive', 'Interpreter', 'latex')
end