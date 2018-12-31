classdef RealDecompositionRep < replab.RealRep
% Describes a direct sum of representations, where all components are all
% pairwise "fully inequivalent"
%
% We say that rho and tau are fully inequivalent if, in the irreducible decomposition
% rho = rho_1 + rho_2 + ... and tau = tau_1 + tau_2 + ..., no rho_i is equivalent to a tau_jx
    properties (SetAccess = protected)
        components;
    end
    methods
        function self = RealDecompositionRep(realRep, components)
            nG = realRep.group.nGenerators;
            images = cell(1, nG);
            imagesInv = cell(1, nG);
            for i = 1:nG
                img = cellfun(@(c) c.images{i}, components, 'UniformOutput', false);
                imgInv = cellfun(@(c) c.imagesInv{i}, components, 'UniformOutput', false);
                images{i} = blkdiag(img{:});
                imagesInv{i} = blkdiag(imgInv{:});
            end
            [U Uinv] = replab.rep.computeU(components);
            if isempty(U)
                parent = [];
            else
                parent = realRep;
            end
            self = self@replab.RealRep(realRep.group, realRep.dimension, images, imagesInv, parent, U, Uinv);
            self.components = components;
        end
        function n = nComponents(self)
            n = length(self.components);
        end
        function c = component(self, r)
            c = self.components{r};
        end
        function s = str(self)
            s = 'Irreducible decomposition with components';
            sep = ' ';
            for i = 1:self.nComponents
                s = sprintf('%s%sI(%d)x%s(%d)', s, sep, self.component(i).multiplicity, self.component(i).divisionAlgebra.shortName, self.component(i).dimension1);
                sep = ' + ';
            end
        end
        function c = centralizerAlgebra(self)
            c = replab.RealDecompositionCentralizerAlgebra(self);
        end
    end
end
