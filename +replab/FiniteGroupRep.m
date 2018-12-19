classdef FiniteGroupRep < replab.Rep
    
    properties (SetAccess = protected)
        images;
    end
    
    methods
        
        function self = FiniteGroupRep(group, images, T)
            assert(length(images) == group.nGenerators);
            self.T = T;
            self.group = group;
            self.images = images;
        end
        
        function d = dimension(self)
            d = self.T.n;
        end
        
        function f = field(self)
            f = self.T.field;
        end

        function rho = image(self, g)
            word = self.group.factorization(g);
            rho = self.T.identity;
            for i = 1:length(word.indices)
                g = self.images{word.indices(i)};
                e = word.exponents(i);
                ge = self.T.composeN(g, e);
                rho = self.T.compose(rho, ge);
            end
        end
        
        function s = str(self)
            s = sprintf('Representation of dimension %d in %s', self.dimension, self.field);
        end
        
        function disp(self)
            disp(self.str);
        end
        
    end
    
    properties (Access = protected)
        centralizerAlgebra_ = [];
        isoDec_ = [];
        blockStructure_ = [];
    end
    
    methods
        
        function I = isoDec(self)
            if isempty(self.isoDec_)
                self.isoDec_ = replab.rep.IsoDec.forAlgebra(self.centralizerAlgebra);
            end
            I = self.isoDec_;
        end
        
        function n = nIsotypicComponents(self)
        % Returns the number of isotypic components in this representation
            n = self.isoDec.nComponents;
        end
        
        function [subrho U dim mul] = isotypicComponent(self, i)
        % Returns the i-th isotypic component
        % Outputs:
        % subrho is the subrepresentation in this isotypic component
        % U      is the change of basis matrix such that subrho = U'*rho*U
        % dim    is the representation dimension
        % mul    is the representation multiplicity 
            U = self.isoDec.compBasis(i);
            subImages = cell(1, length(self.images));
            for i = 1:length(self.images)
                subImages{i} = U'*self.images{i}*U;
            end
            md = size(U, 2);
            subT = replab.GeneralLinearGroup(md, 'R15'); %TODO
            subrho = replab.FiniteGroupRep(self.group, subImages, subT);
            dim = self.isoDec.repDims(i);
            mul = self.isoDec.repMuls(i);
        end
        
        function M = irrepMultiplicities(self)
        end
        
        function m = irrepMultiplicity(self, i)
        end
        
        function D = irrepDimensions(self)
        end
        function d = irrepDimension(self, i)
        end
        
        function [subrho U] = irrep(self, i, j)
        end
        
        function a = centralizerAlgebra(self)
            if isempty(self.centralizerAlgebra_)
                self.centralizerAlgebra_ = replab.rep.Algebra.forRep(self);
            end
            a = self.centralizerAlgebra_;
        end
        
        function p = blockStructure(self)
            if isempty(self.blockStructure_)
                d = self.dimension;
                mask = false(d, d);
                for i = 1:length(self.images)
                    mask = mask | (self.images{i} ~= 0);
                end
                self.blockStructure_ = replab.Partition.connectedComponents(mask);
            end
            p = self.blockStructurex_;
        end
        
        

        
% $$$         function M = sampleGroupElement(self)
% $$$             M = self.image(self.group.sample);
% $$$         end
% $$$ 
% $$$         function M = sampleGroupAlgebra(self, nSamples, nRounds)
% $$$             if nargin < 3
% $$$                 nRounds = 4;
% $$$              end
% $$$              if nargin < 2
% $$$                  nSamples = 5;
% $$$             end                
% $$$             M = self.T.identity;
% $$$             for i = 1:nRounds
% $$$                 M1 = randn * self.image(self.group.sample);
% $$$                 for j = 2:nSamples
% $$$                     M1 = M1 + randn * self.image(self.group.sample);
% $$$                 end
% $$$                 M = M * M1;
% $$$             end
% $$$             M = M / sqrt(nSamples^nRounds);
% $$$         end
        
    end

end
