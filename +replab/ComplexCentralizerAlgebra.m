classdef ComplexCentralizerAlgebra < replab.Domain
% Describes an algebra of matrices that commute with the representation
% of a group
    
    properties (SetAccess = protected)
        complexRep; % representation of which this is the centralizer algebra
        n; % Size of the n x n matrix representation of this algebra
    end
    
    properties (Access = protected)
        transversalImages_ = [];
    end
    
    methods
        
        function self = ComplexCentralizerAlgebra(complexRep)
            assert(isa(complexRep, 'replab.ComplexRep'));
            self.complexRep = complexRep;
            self.n = complexRep.dimension;
        end
        
        function b = eqv(self, X, Y)
        % TODO: replace by something better
            b = replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
        end
        
        function M = sample(self)
            M = sampleGeneric(self); % for the domain sample, use the generic sampling
        end
        
        function M = sampleGeneric(self)
        % Samples a generic matrix from the algebra, with additional guarantees
        %
        % The genericity comes from good separation of eigenvalues
        % TODO: define this genericity formally
            M = replab.rep.sampleComplexMatrix(self.n, self.n);
            M = self.project(M);
        end
        
        function M = sampleGenericSelfAdjoint(self)
        % Samples a generic self-adjoint matrix from the algebra
            M = replab.rep.sampleHermitianMatrix(self.n);
            M = self.project(M);
            M = (M + M')/2;
        end
        
        function T = transversalImages(self)
        % Computes the transveral images, used to speed up the computation
        % of the projection into the algebra
            if isempty(self.transversalImages_)
                T = self.complexRep.group.decomposition.transversals;
                tFun = @(t) cellfun(@(x) self.complexRep.image(x), t, ...
                                    'UniformOutput', false);
                self.transversalImages_ = cellfun(tFun, T, 'UniformOutput', false);
            end
            T = self.transversalImages_;
        end
        
        function M1 = project(self, M)
        % Projects the matrix "M" in this centralizer algebra
        % For centralizer algebras constructed from unitary group representations, this
        % corresponds to an orthogonal projection
            T = self.transversalImages;
            M1 = M;
            for i = length(T):-1:1
                t = T{i};
                S = zeros(size(M1));
                for j = 1:length(t)
                    S = S + t{j}*M1*t{j}';
                end
                M1 = S / length(t);
            end
        end
        
    end
    
end
