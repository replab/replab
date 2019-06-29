classdef Commutant < replab.Domain
% Describes the algebra of n x n matrices that commute with a group representation
%
% Let rep be a representation of a group G. This describes the set
% of matrices X such that rep.image(g) * X = X * rep.image(g)

    properties (SetAccess = protected)
        field; % vector space field, 'R' (real) or 'C' (complex)
        n;     % matrix size
        group; % group
        rep;   % representation
    end
    
    properties (Access = protected)
        parent_; % parent domain: real or complex matrices
    end
    
    methods
        
        function self = Commutant(rep)
            self.rep = rep;
            self.n = rep.dimension;
            self.field = rep.field;
            self.group = rep.group;
            switch self.field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(self.n, self.n);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(self.n, self.n);
              otherwise
                error('Unknown field');
            end
        end
        
        function X = project(self, X)
        % Projects any n x n matrix in the equivariant subspace
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
                gX = self.rep.adjointAction(g, X);
                if i == 1
                    X1 = gX;
                else
                    X1 = X1 + gX;
                end
            end
            X1 = X1/nElements;
        end
        
        function X1 = randomAveraging(self, X, nSamples)
        % Computes the average of X using Monte Carlo sampling
            samples = arrayfun(@(i) self.group.sample, 1:nSamples, 'uniform', 0);
            X1 = self.averageOver(X, samples);
        end
        
        function distances = randomDistances(self, X, nSamples)
        % Samples random group elements, and computes the corresponding violation of commutativity

            distances = [];
            for i = 1:nSamples
                g = self.group.sample;
                gX = self.rep.adjointAction(g, X);
                distances(i) = norm(X - gX, 'fro');
            end
        end

        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d %s commutant matrices', ...
                        self.n, self.n, ...
                        replab.str.realComplex(self.field));
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.project(self.parent_.sample);
        end
        
    end
    
end
