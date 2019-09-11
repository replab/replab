classdef Commutant < replab.Domain
% Describes the algebra of n x n matrices that commute with a group representation
%
% Let rep be a representation of a group G. This describes the set
% of matrices X such that rep.image(g) * X = X * rep.image(g)
%
% Abstract base class

    properties (SetAccess = protected)
        parent; % parent domain: either real/complex matrices
        field;  % vector space field, 'R', 'C'
        n;      % matrix size
        group;  % group
        rep;    % representation
    end
    
    methods % ABSTRACT
        
        function X = project(self, X)
        % Projects any n x n matrix in the invariant subspace
            error('Not implemented');
        end

    end
    
    methods (Access = protected)
        
        function self = Commutant(rep)
        % Constructor; please do not call this from user code, but
        % rather use `rep.commutant` which can select an optimized implementation 
        % depending on the use case.
            self.rep = rep;
            self.n = rep.dimension;
            self.field = rep.field;
            self.group = rep.group;
            switch self.field
              case 'R'
                self.parent = replab.domain.RealMatrices(self.n, self.n);
              case 'C'
                self.parent = replab.domain.ComplexMatrices(self.n, self.n);
              otherwise
                error('Unknown field');
            end
        end
        
    end
    
    methods

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
                        replab.str.field(self.field));
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end
        
        function X = sample(self)
        % Samples a generic matrix from this commutant algebra
            X = self.project(self.parent.sample);
        end
        
        function X = sampleSelfAdjoint(self)
        % Samples a generic self-adjoint matrix from this commutant algebra
        %
        % Enforces that X = X'
        % i.e. real matrices are symmetric, complex matrices are Hermitian
            X = self.parent.sample;
            X = (X + X')/2;
            X = self.project(X);
            X = (X + X')/2;
        end
        
        function X = sampleSkewAdjoint(self)
        % Samples a generic skew-adjoint (or anti-self-adjoint) matrix
        %
        % Enforces that X = -X'
            X = self.parent.sample;
            X = (X - X')/2;
            X = self.project(X);
            X = (X - X')/2;
        end

    end
    
end
