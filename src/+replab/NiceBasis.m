classdef NiceBasis < replab.Str
% Represents a decomposition of a orthogonormal basis into an integer matrix with the same span and a correction factor
%   
% The represented matrix is the product ``U = T * V``.

    properties (SetAccess = protected)
        T % double matrix, can be sparse
        V % integer matrix, can be sparse
    end

    methods

        function self = NiceBasis(T, V)
            self.T = T;
            self.V = V;
        end

        function U = U(self)
            U = self.T * self.V;
        end
        
        function d = childDimension(self)
        % Returns the dimension of the child subrepresentation
        %
        % Returns:
        %   integer: Dimension of the subrepresentation
            d = size(self.V, 1);
        end
        
        function d = parentDimension(self)
        % Returns the dimension of the parent representation
        %
        % Returns:
        %   integer: Dimension of the parent representation
            d = size(self.V, 2);
        end
        
        function b = hasCorrection(self)
        % Returns whether the represented basis has any correction factors
        %
        % Returns:
        %   logical: Whether the correction matrix `T` differs from the identity
            b = isequal(self.T, eye(self.childDimension));
        end
        
        function b = isCorrectionDiagonal(self)
        % Returns whether the represented basis has a diagonal correction matrix
        %
        % Returns:
        %   logical: Whether the correction matrix `T` is diagonal
            b = isdiag(self.T);
        end
        
        function D = normalizationFactors(self)
        % Returns a string representation of the correction coefficients, provided the correction matrix is diagonal
        %
        % Returns:
        %   row cell array of char: String representation of the normalization factor to be appended to a vector
        %                           e.g. the coefficient ``1/2`` is represented by '/2'.
        %
        % Raises:
        %   An error if the correction matrix is not diagonal
            assert(self.isCorrectionDiagonal);
            d = self.childDimension;
            D = cell(1, d);
            for i = 1:d
                c = self.T(i, i);
                inv2 = 1/(c*c);
                if abs(inv2 - round(inv2)) < replab.Settings.doubleEigTol
                    n = round(inv2);
                    if n == 1
                        D{i} = '';
                    elseif sqrt(n) == round(sqrt(n))
                        n1 = sqrt(n);
                        D{i} = sprintf('/%d', n1);
                    else
                        D{i} = sprintf('/sqrt(%d)', n);
                    end
                else
                    D{i} = sprintf('* %e', c);
                end
            end
        end

        function res = mtimes(self, rhs)
        % Returns the decomposition of the composition of two bases
        %
        % We have that the parent representation has a subrepresentation given by the basis
        % expressed by `rhs`, and that subrepresentation has a subrepresentation given by the
        % basis of `self`.
        %  
        % Args:
        %   rhs (replab.NiceBasis): The basis in the middle representation
        %
        % Returns:
        %   replab.NiceBasis: The decomposition of ``self.U * rhs.U``
            if ~rhs.hasCorrection
                res = replab.NiceBasis(self.T, self.V * rhs.V);
            else
                res = replab.NiceBasis.fromV(self.V * rhs.V);
            end
        end
        
    end
    
    methods (Static)

        function niceBasis = fromIntegerBasis(V)
        % Constructs a nice basis from an integer matrix
        %
        % Args:
        %   V (integer matrix): Integer basis
            assert(isreal(V) && ~any(any(isnan(V))));
            VV = V*V';
            if isdiag(VV)
                T = diag(1./sqrt(diag(VV)));
            else
                T = inv(chol(V*V', 'lower'));
            end
            niceBasis = replab.NiceBasis(T, V);
        end
        
        function niceBasisOpt = attemptFromUnitary(U)
        % Attempts integer basis recovery from the given unitary basis
        %
        % Tries to find a rational matrix that has the same row span as the given matrix `U`. It uses the `rat`
        % Matlab function that uses truncated continued fraction expansions. Additionally, the function attemps
        % to perform Gram Schmidt orthogonalization over the integers (but will skip that step if it would make the
        % coefficients to grow too big).
        %
        % The method fails when the recovered approximation has integer coefficients that are too big.
        %
        % Args:
        %   U (double matrix): Unitary/orthonormal (real) basis of dimension dChild x dParent with dChild >= dParent
        %
        % Returns:
        %   replab.NiceBasis or []: Nice decomposition of the basis, or [] if the recovery failed
        %
        % Raises:
        %   An error if the basis is complex
            assert(isreal(U), 'The basis must be real');
            P = U'*U;
            [~, jb] = rref(P);
            U = P(jb, :);
            [num den] = replab.NiceBasis.attemptRecoverRational(U);
            if isempty(num)
                niceBasisOpt = [];
            else
                afterGS = replab.NiceBasis.integerGramSchmidt(num);
                if isempty(afterGS)
                    niceBasisOpt = replab.NiceBasis.fromIntegerBasis(num);
                else
                    niceBasisOpt = replab.NiceBasis.fromIntegerBasis(afterGS);
                end
            end
        end
        
    end

end
