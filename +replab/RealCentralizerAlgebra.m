classdef RealCentralizerAlgebra < replab.Str
    
    properties (SetAccess = protected)
        realRep; % representation of which this is the centralizer algebra
        n; % Size of the n x n matrix representation of this algebra
    end
    
    properties (Access = protected)
        transversalImages_ = [];
    end
    
    methods
        
        function self = RealCentralizerAlgebra(realRep)
            assert(isa(realRep, 'replab.RealRep'));
            self.realRep = realRep;
            self.n = realRep.dimension;
        end
        
        function M = sampleUniformly(self)
            M = replab.rep.sampleRealMatrix(self.n, self.n);
            M = self.project(M);
        end
        
        function M = sampleUniformlySelfAdjoint(self)
            M = replab.rep.sampleSymmetricMatrix(self.n);
            M = self.project(M);
            M = (M + M')/2;
        end
        
        function T = transversalImages(self)
            if isempty(self.transversalImages_)
                T = self.realRep.group.decomposition.transversals;
                tFun = @(t) cellfun(@(x) self.realRep.image(x), t, ...
                                    'UniformOutput', false);
                self.transversalImages_ = cellfun(tFun, T, 'UniformOutput', false);
            end
            T = self.transversalImages_;
        end
        
        function M1 = project(self, M)
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
