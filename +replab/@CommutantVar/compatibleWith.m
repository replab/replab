function okLevel = compatibleWith(X, Y)
% checks whether the block structure of Y is compatible with that of X
%
% the output is:
%  0 if Y is not compatible with X
%  1 if Y has the block structure of X
%  2 if Y has the block structure of X and identical blocks in X
%    correspond to identical blocks in Y

    % Numerical tolerance to decide whether numbers are close to zero in
    % this function.
    epsilon = 1e-10;
    epsilonWarning = 1e-15;

    % We keep track if some rounding is done
    maxOuterEpsilonFound = 0;
    maxEpsilonFound = 0;
    
    % basic tests
    if ~isequal(size(X), size(Y))
        okLevel = 0;
        return;
    end
    
    % We examine each case independently
    okLevel = 2;
    if isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % CommutantVar vs CommutantVar
        
        % We put Y in the block basis of self
        rotatedY = X.U'*Y.fullMatrix*X.U;
        
        % We check the structure of rotatedY
        co = 0;
        for i = 1:X.nComponents
            dim = X.dimensions1(i);
            switch X.types(i)
                case 'R'
                case 'C'
                    dim = dim/2;
                case 'H'
                    dim = dim/4;
                otherwise
                    error('Unknown type');
            end
            
            % First we check the gross structure
            shouldBeZero = rotatedY(co + (1:dim*size(X.blocks{i},1)), co + 1 + dim*size(X.blocks{i},2):end);
            indices = getvariables(shouldBeZero);
            for ind = indices
                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                    okLevel = 0;
                    break;
                end
                maxOuterEpsilonFound = max(maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
            end
            
            % Now we check the intro-block structure
            if okLevel == 2
                % Check full compatibility
                ideal = rotatedY(co + (1:dim:dim*size(X.blocks{i},1)), co + (1:dim:dim*size(X.blocks{i},2)));
                ideal = kron(ideal, eye(dim));
                shouldBeZero = ideal - rotatedY(co + (1:dim*size(X.blocks{i},1)), co + (1:dim*size(X.blocks{i},2)));
                
                % we check if the coefficients are negligible
                indices = getvariables(shouldBeZero);
                for ind = indices
                    if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                        okLevel = 1;
                        maxEpsilonFound = 0;
                        break;
                    end
                    maxEpsilonFound = max(maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                end
            end
            if okLevel == 1
                % Check partial compatibility
                mask = kron(ones(size(X.blocks{i})), eye(dim));
                shouldBeZero = rotatedY(co + (1:dim*size(X.blocks{i},1)), co + (1:dim*size(X.blocks{i},2))).*(1-mask);

                % we check if the coefficients are negligible
                indices = getvariables(shouldBeZero);
                for ind = indices
                    if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                        okLevel = 0;
                        return;
                    end
                    maxEpsilonFound = max(maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                end
            end
            co = co + dim*size(X.blocks{i},1);
        end
    elseif isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar')
        % CommutantVar vs sthg
        
        % We put Y in the block basis of self
        rotatedY = X.U'*Y*X.U;

        % We check the structure of rotatedY
        co = 0;
        for i = 1:X.nComponents
            dim = X.dimensions1(i);
            switch X.types(i)
                case 'R'
                case 'C'
                    dim = dim/2;
                case 'H'
                    dim = dim/4;
                otherwise
                    error('Unknown type');
            end
            
            % First we check the gross structure
            shouldBeZero = rotatedY(co + (1:dim*size(X.blocks{i},1)), co + 1 + dim*size(X.blocks{i},2):end);
            if max(max(abs(shouldBeZero))) > epsilon
                okLevel = 0;
                return;
            end
            maxOuterEpsilonFound = max(maxOuterEpsilonFound, max(max(abs(shouldBeZero))));
            
            % Now we check the intro-block structure
            if okLevel == 2
                % Check full compatibility
                ideal = rotatedY(co + (1:dim:dim*size(X.blocks{i},1)), co + (1:dim:dim*size(X.blocks{i},2)));
                ideal = kron(ideal, eye(dim));
                shouldBeZero = ideal - rotatedY(co + (1:dim*size(X.blocks{i},1)), co + (1:dim*size(X.blocks{i},2)));
                if max(max(abs(shouldBeZero))) > epsilon
                    okLevel = 1;
                    maxEpsilonFound = 0;
                end
                maxEpsilonFound = max(maxEpsilonFound, max(max(abs(shouldBeZero))));
            end
            if okLevel == 1
                % Check partial compatibility
                mask = kron(ones(size(X.blocks{i})), eye(dim));
                shouldBeZero = rotatedY(co + (1:dim*size(X.blocks{i},1)), co + (1:dim*size(X.blocks{i},2))).*(1-mask);
                if max(max(abs(shouldBeZero))) > epsilon
                    okLevel = 0;
                    return;
                end
                maxEpsilonFound = max(maxEpsilonFound, max(max(abs(shouldBeZero))));
            end
            co = co + dim*size(X.blocks{i},1);
        end
    elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % sthg vs CommutantVar
        okLevel = Y.compatibleWith(X);
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
    
    
    % We produce a warning if some small but not too small coefficients
    % have been neglected
    if max(maxOuterEpsilonFound, maxEpsilonFound) > epsilonWarning
        warning(['Block structure mismatch by ', num2str(maxEpsilonFound), ' ignored.']);
    end
end
