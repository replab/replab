function Z = rdivide(X, Y)
    % element-wise right division
    
    % Check that dimensions are compatible
    if ~isscalar(Y)
        error('Use fullMatrix for non-scalar multiplications.');
    end
    if ~isa(X, 'replab.CommutantVar') || isa(Y, 'replab.CommutantVar')
        error('Numerator should be a CommutantVar and denominator should not.');
    end
    
    Z = X.copy;
    for i = 1:Z.nComponents
        Z.blocks{i} = Z.blocks{i}./Y;
    end
end
