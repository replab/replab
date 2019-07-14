classdef Equivariant < replab.Domain
% Describes a vector space of group-equivariant matrices
%
% Let repR and repC be two representations of the same group G.
%
% This describes the set of matrices X such that repR.image(g) * X = X * repC.image(g)
%
% See Proposition 4 of 
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
    
    properties (SetAccess = protected)
        field; % vector space field, 'R' (real) or 'C' (complex)
        nR;    % row size
        nC;    % column size
        group; % group
        repR;  % representation of row space
        repC;  % representation of column space
    end
    
    properties (Access = protected)
        parent_; % parent domain: real or complex matrices
    end
    
    methods

        function self = Equivariant(repR, repC)
        % Constructor; please do not call this from user code, but
        % rather use `replab.rep.equivariant(repR, repC)`, 
        % or `repR.equivariant(repC)`, which can eventually
        % select an optimization implementation depending on the
        % use case.
            self.repR = repR;
            self.nR = repR.dimension;
            self.repC = repC;
            self.nC = repC.dimension;
            assert(isequal(repR.field, repC.field), ...
                   'Both representations must have be defined on the same field');
            self.field = repR.field;
            assert(isequal(repR.group, repC.group), ...
                   'Both representations must be defined on the same group');
            self.group = repR.group;
            switch self.field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(self.nR, self.nC);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(self.nR, self.nC);
              case 'H'
                self.parent_ = replab.domain.QuaternionMatrices(self.nR, self.nC);                
              otherwise
                error('Unknown field');
            end
        end
        
        function X = project(self, X)
        % Projects any nR x nC matrix in the equivariant subspace
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
                gX = self.repC.action(g, self.repR.action(g, X)')';
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
        
        function distances = sampleDistances(self, X, nSamples)
        % Samples random group elements, and computes the corresponding violation of equivariance
            distances = [];
            for i = 1:nSamples
                g = self.group.sample;
                gX = self.repC.action(g, self.repR.action(g, X)')';
                distances(i) = norm(X - gX, 'fro');
            end
        end

        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d %s equivariant matrices', ...
                        self.nR, self.nC, ...
                        replab.str.field(self.field));
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
