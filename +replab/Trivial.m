classdef Trivial < replab.Domain
% Describes the trivial subrepresentation of a group representation
    properties (SetAccess = protected)
        field;
        n;
        group;
        rep;
    end
        
    properties (Access = protected)
        parent_; % parent domain: real or complex vectors
    end
    
    methods
        
        function self = Trivial(rep)
            self.rep = rep;
            self.field = rep.field;
            self.n = rep.dimension;
            self.group = rep.group;
            switch self.field
              case 'R'
                self.parent_ = replab.domain.RealVectors(self.n);
              case 'C'
                self.parent_ = replab.domain.ComplexVectors(self.n);
              otherwise
                error('Unknown field');
            end
        end

        function b = isInvariant(self, x)
            b = self.parent_.eqv(self.project(x), x);
        end
        
        function X = project(self, X)
        % Projects vectors in the invariant subspace
        %
        % One can provide n vectors to that function, stored as a
        % d x n matrix
            assert(isa(self.group, 'replab.FiniteGroup'));
            T = self.group.decomposition.transversals;
            for i = length(T):-1:1
                X = self.averageOver(X, T{i});
            end
        end

        function X1 = averageOver(self, X, elements)
        % Averages X over the given group elements (row cell vector)
            X1 = [];
            nElements = length(elements);
            for i = 1:nElements
                g = elements{i};
                gX = self.rep.action(g, X);
                if i == 1
                    X1 = gX;
                else
                    X1 = X1 + gX;
                end
            end
            X1 = X1/nElements;
        end

    end

end
