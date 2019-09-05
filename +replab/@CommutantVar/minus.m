function Z = minus(X,Y)
    % substraction operator
    
    size1 = size(X);
    size2 = size(Y);

    % Check that dimensions are compatible
    if ~isequal(size1, size2)
        error('Incompatible size for matrix substraction');
    end

    % We verify that both variables have compatible structures
    compatLevel = X.compatibleWith(Y);
    if compatLevel == 0
        error('Block structure of both matrices don''t match. Consider using fullMatrix.');
    end

	% We examine each case independently
    if isa(X, 'replab.CommutantVar') && isa(Z, 'replab.CommutantVar')
        % CommutantVar - CommutantVar
    elseif isa(X, 'replab.CommutantVar') && ~isa(Z, 'replab.CommutantVar')
        % CommutantVar - sthg
    elseif ~isa(X, 'replab.CommutantVar') && isa(Z, 'replab.CommutantVar')
        % sthg - CommutantVar
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
    
    Z = X.copy;
    for i = 1:Y.nComponents
        Y.blocks{i} = -Y.blocks{i};
    end
end
