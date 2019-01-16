classdef PermutationGroup < replab.FiniteGroup
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
        
        function rho = naturalRepresentation(self)
            rho = self.permutationRepresentation(self.domainSize, self.generators);
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
