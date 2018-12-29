classdef SignedPermutationGroup < replab.FiniteGroup
% Represents a signed permutation group

    properties (SetAccess = protected)
        domainSize; % d when this group acts on {-d..-1, 1..d}
    end
    
    properties (Access = protected)
        abs_ = [];
    end
    
    methods

        function self = SignedPermutationGroup(domainSize)
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
        end
        
        % Domain
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        % Semigroup/Monoid/Group
        
        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end
        
        function y = inverse(self, x)
            n = self.domainSize;
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end
        
        % Own methods
                
        function A = naturalAction(self)
        % Returns the action of elementso f this group on
        % its domain {-d..-1 1..d} where d is self.domainSize
            A = replab.perm.SignedPermutationNaturalAction(self);
        end
        
        function A = vectorAction(self)
        % Returns the action of elements of this group on
        % (self.domainSize)-dimensional vectors by permuting their
        % coefficients and flipping their signs
            A = replab.perm.SignedPermutationVectorAction(self);
        end

        function A = matrixAction(self)
        % Returns the action of elements of this group on d x d matrices
        % where d = self.domainSize, by simultaneous permutations of their
        % rows and columns and flipping their signs
            A = replab.perm.SignedPermutationMatrixAction(self);
        end
        
        function rho = naturalRepresentation(self)
        % Natural representation of signed permutations on integer -d..-1, 1..d
            rho = self.signedPermutationRepresentation(self.domainSize, self.generators);
        end
        
        function G = abs(self)
            if isempty(self.abs_)
                % maps generators to their "absolute value", collapsing
                % i and -i in the domain, and removes generators that map
                % to the identity
                absGenerators = cell(1, 0);
                for i = 1:self.nGenerators
                    ag = abs(self.generators{i});
                    if ~self.isIdentity(ag)
                        absGenerators{end+1} = ag;
                    end
                end
                self.abs_ = replab.Permutations(self.domainSize).subgroup(absGenerators);
            end
            G = self.abs_;
        end
        
    end

    methods (Static)
        
        function perm = toPermutation(signedPerm)
        % Returns the permutation corresponding to the given signed permutation
        % where the permutation acts on the list [-d,..,-1,1,..,d]
            n = length(signedPerm);
            perm = zeros(1, 2*n);
            for i = 1:length(signedPerm)
                im = n - i + 1; % position of -i in the list
                ip = n + i; % position of i in the list
                j = signedPerm(i);
                jm = n - abs(j) + 1;
                jp = n + abs(j);
                if j > 0
                    perm(im) = jm;
                    perm(ip) = jp;
                else
                    perm(im) = jp;
                    perm(ip) = jm;
                end
            end
        end
        
        function signedPerm = fromPermutation(perm)
        % Returns the signed permutation corresponding to the given permutation encoding
        % see SignedPermutations.toPermutation
            n2 = length(perm);
            if mod(n2, 2) ~= 0
                error('Not an image of a signed permutation');
            end
            n = n2/2; % domain size of the signed permutation
            perm(perm <= n) = -(n - perm(perm <= n) + 1);
            perm(perm > n) = perm(perm > n) - n;
            mperm = -fliplr(perm(1:n));
            pperm = perm(n+1:n2);
            assert(isequal(mperm, pperm), 'Not an image of a signed permutation');
            signedPerm = pperm;
        end

        function mat = toMatrix(signedPerm)
        % Returns the signed permutation matrix corresponding to the given
        % signed permutation such that matrix multiplication is
        % compatible with composition of permutations, i.e. 
        % S.toMatrix(S.compose(x, y)) = 
        % S.toMatrix(x) * S.toMatrix(y)
        % where S = SignedPermutations(domainSize)
            n = length(signedPerm);
            mat = sparse(abs(signedPerm), 1:n, sign(signedPerm), n, n);
        end
        
        function signedPerm = fromMatrix(mat)
        % Returns the signed permutation corresponding to the given matrix representation
        % or throws an error
            if isequal(size(mat), [0 0])
                signedPerm = zeros(1, 0);
                return
            end
            signedPerm = [];
            n = size(mat, 1);
            [I J S] = find(mat);
            if length(I) ~= n
                error('Not a monomial matrix');
            end
            I = I';
            J = J';
            S = S';
            sI = sort(I);
            [sJ IJ] = sort(J);
            if ~isequal(sI, 1:n) || ~isequal(sJ, 1:n)
                error('Not a monomial matrix');
            end
            if ~isequal(abs(S), ones(1, n))
                error('Monomial matrix with entries other than +1, -1');
            end
            signedPerm = I.*S;
            signedPerm = signedPerm(IJ);
        end
    
    end
    
end
