function FLAG = isequal_approx(x1, x2, tol)
    % Test that to values are approximately equal given a tolerance
    FLAG = all(abs(real(x1(:) - x2(:))) < tol) & all(abs(imag(x1(:) - x2(:))) < tol);
end