classdef IsotypicDecomposition < replab.Str
% Isotypic decomposition of a representation
    
    properties (SetAccess = protected)
        rep; % Representation being decomposed
        B;   % Block diagonal matrix algebra
        nComponents;   % Number of isotypic components
    end
    
    properties
        isoDec; % Wrapper of IsoDec object during migration
    end

    methods
        
        function self = IsotypicDecomposition(rep, isoDec)
            self.rep = rep;
            self.isoDec = isoDec;
            blockSizes = isoDec.compDims;
            blocks = arrayfun(@(n) replab.Matrices(n, n, false, false), blockSizes, 'UniformOutput', false);
            self.B = replab.BlockDiagonalMatrices(blocks);
            self.nComponents = isoDec.nComponents;
        end
        
        function s = str(self)
            sz = num2str(self.isoDec.compDims);
            s = sprintf('Isotypic decomposition with %d components of sizes %s', self.nComponents, sz);
        end
        
        function I = withTrivialFirst(self)
            nC = self.nComponents;
            nG = length(self.rep.images);
            dev = zeros(nG, nC);
            for i = 1:nC
                d = self.isoDec.compDims(i);
                U = self.isoDec.compBasis(i);
                for j = 1:nG
                    dev(i, j) = norm(U'*self.rep.images{j}*U - eye(d));
                end
            end
            [~, trivial] = min(sum(dev.^2, 2));
            rest = setdiff(1:nC, trivial);
            p = [trivial rest];
            I = replab.IsotypicDecomposition(self.rep, self.isoDec.permuteComponents(p));
        end
        
        function I1 = refined(self)
            I1 = replab.IsotypicDecomposition(self.rep, self.isoDec.refined);
        end
        
        function U = adaptedBasis(self)
            U = self.isoDec.U;
        end
        
        function [subRep U] = component(self, i)
        % Returns the i-th isotypic component
        % Outputs:
        % subRep is the subrepresentation in this isotypic component
        % U      is the change of basis matrix such that subRep = U'*rep*U
            U = self.isoDec.compBasis(i);
            subImages = cellfun(@(X) U'*X*U, self.rep.images, 'UniformOutput', false);
            subd = size(U, 2);
            subT = replab.GeneralLinearGroup(subd, false);
            isUnitary = true; %TODO
            subRep = replab.FiniteGroupRep(self.rep.group, subImages, isUnitary, subT);
        end
        
        function blocks = blocksOf(self, M)
        % Let M be a matrix block-diagonal in the isotypic basis,
        % i.e. M \in  V1 x V1' + V2 x V2' + ... where Vi are isotypic subspaces
        %
        % This returns a decomposition of M into blocks of size
        % dim(V1)xdim(V1), dim(V2)xdim(V2), etc..
            blocks = cell(1, self.nComponents);
            for i = 1:self.nComponents
                U = self.isoDec.compBasis(i);
                blocks{i} = U' * M * U;
            end
        end
        
    end
    
    methods (Static)
        
        function I = ofRep(rep)
            isoDec = replab.rep.IsoDec.fromAlgebra(rep.centralizerAlgebra);
            I = replab.IsotypicDecomposition(rep, isoDec).withTrivialFirst;
        end
    end            
    
end
