function unit_test_linear_system(u, v, w, u_dilatational, v_dilatational, w_dilatational, dx, dy, dz)
    % Unit Test

    % Compute divergence of the velocity 
    [grad_u, grad_v, grad_w] = gradient_velocity_yz_periodic(u, v, w, dx, dy, dz);

    div_1 = grad_u + grad_v + grad_w;
    
    % Compute divergence of the velocity after solving div(grad(phi)) = div(u, v, w)
    [du_dilatational, dv_dilatational, dw_dilatational] = gradient_velocity_yz_periodic(u_dilatational, v_dilatational, w_dilatational, dx, dy, dz);
    [du_solenoidal, dv_solenoidal, dw_solenoidal] = gradient_velocity_yz_periodic(u - u_dilatational, v - v_dilatational, w - w_dilatational, dx, dy, dz);

    div_2 = du_dilatational + dv_dilatational + dw_dilatational + ...
            du_solenoidal + dv_solenoidal + dw_solenoidal;
    
    % Check
    FLAG = isequal_approx(div_1, div_2, 1e-10);
    
    if FLAG
        fprintf('TEST: OK!\n');
    else
        fprintf('TEST: NOT PASSED!\n');
    end

end