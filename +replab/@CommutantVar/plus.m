function Z = plus(X,Y)
    % addition operator
    
    size1 = size(X);
    size2 = size(Y);

    if (~isequal(size1, size2)) && (prod(size1) ~= 1) && (prod(size2) ~= 1) && (prod(size1)+prod(size2) > 0)
        error('Incompatible size for matrix addition');
    end

    % We examine each case independently
    if isa(Y, 'replab.CommutantVar') && isa(Z, 'replab.CommutantVar')
        % CommutantVar + CommutantVar
    elseif isa(Y, 'replab.CommutantVar') && ~isa(Z, 'replab.CommutantVar')
        % CommutantVar + sthg
    elseif ~isa(Y, 'replab.CommutantVar') && isa(Z, 'replab.CommutantVar')
        % sthg + CommutantVar
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
    
    Z = X.copy;
    for i = 1:Y.nComponents
        Y.blocks{i} = -Y.blocks{i};
    end
end
