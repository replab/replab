classdef RealIrreducibleCentralizerAlgebra < replab.RealCentralizerAlgebra
    properties
        nBlocks;
        blockSize;
    end
    methods
        function self = RealIrreducibleCentralizerAlgebra(realIrrep)
            assert(isa(realIrrep, 'replab.RealIrrep'));
            self = self@replab.RealCentralizerAlgebra(realIrrep);
            self.nBlocks = realIrrep.dimension1 / realIrrep.divisionAlgebra.m;
            self.blockSize = realIrrep.multiplicity * realIrrep.divisionAlgebra.m;
        end
        
        function B = block(self, M, c)
        % Returns the c-th block of M, or, if c = 0 or not given, the average over all blocks
            if nargin < 3 || c == 0
                B = self.block(M, 1);
                for c = 2:self.nBlocks
                    B = B + self.block(M, c);
                end
                B = B / self.nBlocks;
                return
            end
            divisionAlgebra = self.realRep.divisionAlgebra;
            % Representation dimension
            d = self.realRep.dimension1;
            % Representation multiplicity
            m = self.realRep.multiplicity;
            % Division algebra dimension
            ad = divisionAlgebra.m;
            % Representation dimension in splitting field
            bs = self.blockSize;
            % Find indices of the c-th block
            ind = bsxfun(@plus, 0:d:d*(m-1), (1:ad)');
            ind1 = ind(:)' + (c-1)*ad;
            B = divisionAlgebra.projectMatrix(M(ind1, ind1));
        end

        function M1 = project(self, M)
            B = self.block(M, 0);
            nB = self.nBlocks;
            bs = self.blockSize;
            B = kron(B, eye(nB));
            ad = self.realRep.divisionAlgebra.m;
            B = reshape(B, [nB ad bs/ad nB ad bs/ad]);
            B = permute(B, [2 1 3 5 4 6]);
            M1 = reshape(B, [nB*bs nB*bs]);
        end
        
        function B = blockOfParentElement(self, M, c)
            if nargin < 3 || c == 0
                B = self.blockOfParentElement(M, 1);
                for c = 2:self.nBlocks
                    B = B + self.blockOfParentElement(M, c);
                end
                B = B / self.nBlocks;
                return
            end
            divisionAlgebra = self.realRep.divisionAlgebra;
            % Representation dimension
            d = self.realRep.dimension1;
            % Representation multiplicity
            m = self.realRep.multiplicity;
            % Division algebra dimension
            ad = divisionAlgebra.m;
            % Representation dimension in splitting field
            bs = self.blockSize;
            % Find indices of the c-th block
            ind = bsxfun(@plus, 0:d:d*(m-1), (1:ad)');
            ind1 = ind(:)' + (c-1)*ad;
            U = self.realRep.U(:, ind1);
            Uinv = self.realRep.Uinv(ind1, :);
            B = Uinv * M * U;
            B = divisionAlgebra.projectMatrix(B);
        end            
        
        function blocks = blocksOfParentElement(self, M)
            blocks = {self.blockOfParentElement(M, 0)};
        end
            
        %        function sampleGeneric(self)
        %        end
        
        %        function sampleGenericSelfAdjoint(self)
        %        end
    end
end
