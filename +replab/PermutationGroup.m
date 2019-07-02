classdef PermutationGroup < replab.Group
% Represents a permutation group

    properties (SetAccess = protected)
        domainSize; % d when this group acts on {1, ..., d}
    end
    
    properties (Access = protected)
        orbits_ = [];
    end
    
    methods
                
        function self = PermutationGroup(domainSize)
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
        end
        
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
        
        function rho = indexPermutationRep(self, localDimension)
        % Representation that permutes the indices of a tensor
        %
        % It acts on the tensor space R^ld x R^ld ... (domainSize times)
        % by permuting the indices. The representation returned is real.
            rho = replab.rep1.IndexPermutationRep(self, localDimension);
        end
        
        function rho = naturalRepresentation(self)
            rho = self.permutationRepresentation(self.domainSize, self.generators);
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
            U = replab.rep1.standardBasis(self.domainSize);
            rho = self.naturalRep.subRep(U(:, 2:end));
        end
        
        function o = orbits(self)
            if isempty(self.orbits_)
                G = zeros(self.nGenerators, self.domainSize);
                for i = 1:self.nGenerators
                    G(i, :) = self.generators{i};
                end
                self.orbits_ = replab.Partition.permutationsOrbits(G);
            end
            o = self.orbits_;
        end
        
    end
        
end
