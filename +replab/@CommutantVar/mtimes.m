function Z = mtimes(X, Y)
    % matrix multiplication
    
    % Check that dimensions are compatible
    if ~isscalar(X) && ~isscalar(Y)
        error('Incompatible size for multiplication, use fullMatrix for non-scalar multiplications.');
    end
    
    Z = X.*Y;
end
