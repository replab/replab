function Z = mldivide(X, Y)
    % matrix left division

    % Only scalar division is supported
    if ~isscalar(X)
        error('Use fullMatrix for non-scalar division.');
    end
    
    Z = Y./X;
end
