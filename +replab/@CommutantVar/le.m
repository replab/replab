function F = le(X,Y)
    % greater or equal constraint

    % We only support some simple cases for now
    if isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % CommutantVar <= CommutantVar
        error('Comparison between two CommutantVar not yet supported');
    elseif isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar')
        % CommutantVar <= sthg
        if isscalar(Y)
            F = [X.blocks{1} <= Y];
            for i = 2:self.nComponents
                F = [F, self.blocks{i} <= Y];
            end
        else
            % We need to check whether X and Y have the same block
            % structure... otherwise, request using fullMatrix mode.
            error('Comparison with a matrix not yet supported')
        end
    elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % sthg <= CommutantVar
        F = ge(Y,X);
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
end
