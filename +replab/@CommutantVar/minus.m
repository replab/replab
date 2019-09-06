function Z = minus(X,Y)
% substraction operator
    
    % Numerical tolerance to decide whether numbers are close to zero in
    % this function.
    epsilon = 1e-10;

    % We keep track of the encountered non-hermiticities
    maxNonHermiticity = 0;
    epsilonWarning = 1e-14;
    
    size1 = size(X);
    size2 = size(Y);

    % Check that dimensions are compatible
    if ~isequal(size1, size2)
        error('Incompatible size for matrix substraction');
    end

	% We examine each case independently
    if isa(X, 'replab.CommutantVar')
        % CommutantVar - sthg
        
        % We verify that both variables have fully compatible structures
        compatLevel = X.compatibleWith(Y);
        if compatLevel ~= 2
            error('Block structure of both matrices don''t match. Consider using fullMatrix.');
        end

        % We express Y in the block basis of X
        if isa(Y, 'replab.CommutantVar')
            rotatedY = X.U'*Y.fullMatrix*X.U;
        else
            rotatedY = X.U'*Y*X.U;
        end

        % The block structure matches fully, we procede to perform the
        % substraction on each block
        Z = copy(X);
        co = 0;
        for i = 1:Z.nComponents
            dim = Z.dimensions1(i);
            switch Z.types(i)
                case 'R'
                case 'C'
                    dim = dim/2;
                case 'H'
                    dim = dim/4;
                otherwise
                    error('Unknown type');
            end
            constantBlock = rotatedY(co + (1:dim:dim*size(X.blocks{i},1)), co + (1:dim:dim*size(X.blocks{i},1)));
            % We make sure that the constant block is hermitian
            shouldBeZero = constantBlock - constantBlock';
            if isa(shouldBeZero, 'sdpvar')
                indices = getvariables(shouldBeZero);
                for ind = indices
                    if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                        error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                    end
                    maxNonHermiticity = max(maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                end
            else
                if max(max(abs(shouldBeZero))) > epsilon
                    error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                end
                maxNonHermiticity = max(maxNonHermiticity, max(max(abs(shouldBeZero))));
            end
            if ~ishermitian(constantBlock)
                % We force exact hermiticity
                constantBlock = (constantBlock + constantBlock')/2;
            end
            Z.blocks{i} = Z.blocks{i} - constantBlock;
            co = co + dim*size(X.blocks{i},1);
        end
    elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
        % sthg - CommutantVar
        Z = -(Y-X);
    else
        error('Neither of the two arguments is of type replab.CommutantVar');
    end
    
    
    % We produce a warning if some small but not too small coefficients
    % have been neglected
    if maxNonHermiticity > epsilonWarning
        warning(['Non-hermiticity of order ', num2str(maxEpsilonFound), ' was corrected.']);
    end
end
