function Z = times(X, Y)
    % element-wise multiplication
    
    % Check that dimensions are compatible
    if ~isscalar(X) && ~isscalar(Y)
        error('Use fullMatrix for non-scalar multiplications.');
    end
    
	% We examine each case independently
    if isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % CommutantVar .* CommutantVar
        error('Use fullMatrix instead');
    elseif isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar')
        % CommutantVar .* sthg
        Z = X.copy;
        for i = 1:Z.nComponents
            Z.blocks{i} = Y.*Z.blocks{i};
        end
    elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % sthg .* CommutantVar
        Z = Y.*X;
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
end
