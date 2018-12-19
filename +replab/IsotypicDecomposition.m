classdef IsotypicDecomposition < replab.Str
% Isotypic decomposition of a representation
    
    properties (SetAccess = protected)
        rep; % Representation being decomposed
        B;   % Block diagonal matrix algebra
        n;   % Number of isotypic components
    end
    
    properties (Access = protected)
        isoDec; % Wrapper of IsoDec object during migration
    end

    methods
        
        function s = str(self)
            sz = num2str(self.isoDec.compDims);
            s = sprintf('Isotypic decomposition with %d components of sizes %s', self.n, sz);
        end
        
        function self = IsotypicDecomposition(rep, isoDec)
            self.rep = rep;
            self.isoDec = isoDec;
            blockSizes = isoDec.compDims;
            blocks = arrayfun(@(n) replab.Matrices(n, n, false, false), blockSizes, 'UniformOutput', false);
            self.B = replab.BlockDiagonalMatrices(blocks);
            self.n = isoDec.nComponents;
        end
        
        function [subRep U] = component(self, i)
        % Returns the i-th isotypic component
        % Outputs:
        % subRep is the subrepresentation in this isotypic component
        % U      is the change of basis matrix such that subrho = U'*rho*U
            U = self.isoDec.compBasis(i);
            subImages = cellfun(@(X) U'*X*U, self.rep.images, 'UniformOutput', false);
            subd = size(U, 2);
            subT = replab.GeneralLinearGroup(subd, false);
            isUnitary = true;
            subRep = replab.FiniteGroupRep(self.rep.group, subImages, isUnitary, subT);
        end
        
        function blocks = blocksOf(self, M)
        % Let M be a matrix block-diagonal in the isotypic basis,
        % i.e. M \in  V1 x V1' + V2 x V2' + ... where Vi are isotypic subspaces
        %
        % This returns a decomposition of M into blocks of size
        % dim(V1)xdim(V1), dim(V2)xdim(V2), etc..
            blocks = cell(1, self.n);
            for i = 1:self.n
                U = self.isoDec.compBasis(i);
                blocks{i} = U' * M * U;
            end
        end
        
    end
    
    methods (Static)
        
        function I = ofRep(rep)
            isoDec = replab.rep.IsoDec.fromAlgebra(rep.centralizerAlgebra, rep.fibers);
            I = replab.IsotypicDecomposition(rep, isoDec);
        end
    end            
    
end
