function Z = mrdivide(X, Y)
    % matrix right division

    % Only scalar division is supported
    if ~isscalar(Y)
        error('Use fullMatrix for non-scalar division.');
    end
    
    Z = X./Y;
end
