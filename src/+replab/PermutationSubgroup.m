classdef PermutationSubgroup < replab.PermutationGroup & replab.NiceFiniteSubgroup
% Represents a permutation group

    properties (SetAccess = protected)
        domainSize; % d when this group acts on {1, ..., d}
    end
    
    methods
        
        % Domain
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        % Semigroup/Monoid/Group
        
        function z = compose(self, x, y)
            z = x(y);
        end
        
        function y = inverse(self, x)
            n = self.domainSize;
            y = zeros(1, n);
            y(x) = 1:n;
        end
        
        function A = naturalAction(self)
        % Returns the action of elements of this group on {1..self.domainSize}
            A = replab.perm.PermutationNaturalAction(self);
        end
        
        function A = vectorAction(self)
        % Returns the action of elements of this group on
        % (self.domainSize)-dimensional vectors by permuting their coefficients
            A = replab.perm.PermutationVectorAction(self);
        end

        function A = matrixAction(self)
        % Returns the action of elements of this group on d x d matrices
        % where d = self.domainSize, by simultaneous permutations of their
        % rows and columns
            A = replab.perm.PermutationMatrixAction(self);
        end
        
        function perm = indexRelabelingPermutation(self, g, localDimension)
        % Describes the permutation action of this group on tensor coefficients
            n = self.domainSize;
            dims = localDimension * ones(1, n);
            perm = permute(reshape(1:prod(dims), dims), fliplr(n +  1 - g));
            perm = perm(:)';
        end
        
        function phi = indexRelabelingMorphism(self, localDimension)
        % Describes the permutation action of this group on tensor coefficients
        %
        % The tensor coefficients correspond to R^ld x R^ld ... (domainSize times)
            phi = @(g) self.indexRelabelingPermutation(g, localDimension);
        end
        
        function rho = indexRelabelingRep(self, localDimension)
        % Representation that permutes the indices of a tensor
        %
        % It acts on the tensor space R^ld x R^ld ... (domainSize times)
        % by permuting the indices. The representation returned is real.
            rho = replab.rep.IndexRelabelingRep(self, localDimension);
        end
        
        function rho = naturalRep(self)
        % Returns the natural permutation representation of this permutation group
            rho = self.permutationRep(self.domainSize, self.generators);
        end
        
        function rho = standardRep(self)
        % Returns the standard representation of this permutation group
        %
        % It corresponds to the representation orthogonal to the
        % trivial representation with basis [1, 1, ..., 1]'/sqrt(d)
            U = replab.rep.standardBasis(self.domainSize);
            U = U(2:end, :);
            rho = self.naturalRep.subRep(U);
        end
        
        function o = orbits(self)
            G = zeros(self.nGenerators, self.domainSize);
            for i = 1:self.nGenerators
                G(i, :) = self.generators{i};
            end
            o = replab.Partition.permutationsOrbits(G);
        end
        
    end
        
end
