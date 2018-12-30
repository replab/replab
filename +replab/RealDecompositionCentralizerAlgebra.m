classdef RealDecompositionCentralizerAlgebra < replab.RealCentralizerAlgebra
% Direct sum of centralizer algebras, coming from inequivalent representations; thus
% "off-diagonal blocks" are zero, by Schur's lemma
    properties (SetAccess = protected)
        components;
    end
    methods
        
        function self = RealDecompositionCentralizerAlgebra(realDecompositionRep)
            assert(isa(realDecompositionRep, 'replab.RealDecompositionRep'));
            self = self@replab.RealCentralizerAlgebra(realDecompositionRep);
            self.components = cellfun(@(x) x.centralizerAlgebra, realDecompositionRep.components, 'UniformOutput', false);
        end
        
        function blocks = blocksOfParentElement(self, M)
            blocks = cell(1, 0);
            for i = 1:length(self.components)
                blocks = horzcat(blocks, self.components{i}.blocksOfParentElement(M));
            end
        end
    end
end
