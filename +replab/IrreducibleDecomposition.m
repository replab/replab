classdef IrreducibleDecomposition < replab.Str
% Irreducible decomposition of a representation

    properties (SetAccess = protected)
        rep; % Representation being decomposed
        U;   % Change of basis matrix
        nComponents;    % Number of isotypic components
        
        % The next properties are given as integer 1 x nComponents vectors,
        % with one number for each irreducible representation/component
        
        dimensions;                % Representation dimension in the current field
        splittingFieldDimensions;  % Representation dimension in a splitting field
                                   % = number of corresponding blocks in the centralizer
        multiplicities;            % Representation multiplicity in the current field
        centralizerBlockDimensions;% Dimension of unique blocks of centralizer elements
                                   % in the current field
        divisionAlgebraDimensions; % Dimension of the division algebra over the
                                   % current field present in the centralizer block

        divisionAlgebras; % Definition of division algebras, as 1 x nComponents cell array

        
        centralizerAlgebraBlocks; % Block diagonal matrix algebra
        groupAlgebraBlocks;       % Block diagonal matrix algebra
    end

    methods
        
        function self = IrreducibleDecomposition(rep, U, dimensions, multiplicities, divisionAlgebras)
            self.rep = rep;
            self.U = U;
            self.nComponents = length(dimensions);
            self.dimensions = dimensions;
            self.multiplicities = multiplicities;
            self.divisionAlgebraDimensions = cellfun(@(A) A.d, divisionAlgebras);
            self.splittingFieldDimensions = self.dimensions ./ self.divisionAlgebraDimensions;
            self.centralizerBlockDimensions = self.multiplicities .* self.divisionAlgebraDimensions;
            blockFun = @(n) replab.Matrices(n, n, false, false);
            self.centralizerAlgebraBlocks = arrayfun(blockFun, self.centralizerBlockDimensions, 'UniformOutput', false);
            self.groupAlgebraBlocks = arrayfun(blockFun, self.dimensions, 'UniformOutput', false);
            self.divisionAlgebras = divisionAlgebras;
        end
        
        function U = adaptedBasis(self)
        % Returns the symmetry adapted change of basis matrix
            U = self.U;
        end
        
        function s = str(self)
            s = 'Irreducible decomposition with components';
            sep = ' ';
            for i = 1:self.nComponents
                s = sprintf('%s%sI(%d)x%s(%d)', s, sep, self.multiplicities(i), self.divisionAlgebras{i}.shortName, self.dimensions(i));
                sep = ' + ';
            end
        end
        
        function [subRep U] = representation(self, r, c)
        % Returns the c-th copy of the r-th representation, where
        % c = 1...self.m(r)
        %
        % If c is omitted or = 0, returns the average of all copies
        %
        % Outputs:
        % subRep is the irreducible subrepresentation
        % U      is the change of basis matrix such that subRep = U'*rep*U
            if nargin < 3 || c == 0
                
                assert(self.multiplicities(r) == 1, 'The copy index can be omitted only for multiplicity 1');
                c = 1;
            end
            subImages = cellfun(@(X) self.blockOfGroupAlgebraMatrix(X, r), self.rep.images, 'UniformOutput', false);
            subT = replab.GeneralLinearGroup(self.dimensions(r), false);
            isUnitary = true; %TODO
            subRep = replab.FiniteGroupRep(self.rep.group, subImages, isUnitary, subT);
        end
        
        function block = blockOfCentralizer(self, M, r, c)
        % A matrix M that commutes with the represenation is split into
        % blocks, with the number of blocks equal to the dimension of the
        % representation
        %
        % This function returns the c-th block of the r-th irreducible representation,
        % where i = 1...splittingFieldDimensions(r)
        %
        % If c is omitted or =0, it returns the average over all copies
            if nargin < 4 || c == 0
                block = self.blockOfCentralizer(M, r, 1);
                for c = 2:self.splittingFieldDimensions(r)
                    block = block + self.blockOfCentralizer(M, r, c);
                end
                block = block / self.splittingFieldDimensions(r);
                return
            end
            sizes = self.multiplicities.*self.dimensions;
            % Size of previous components to skip
            shift = sum(sizes(1:r-1));
            % Representation dimension
            d = self.dimensions(r); 
            % Representation multiplicity
            m = self.multiplicities(r); 
            % Division algebra dimension
            ad = self.divisionAlgebraDimensions(r);
            % Representation dimension in splitting field = number of blocks
            sd = self.splittingFieldDimensions(r); 
            cbd = self.centralizerBlockDimensions(r);
            block = zeros(cbd, cbd);
            % Find indices of the c-th block
            ind = bsxfun(@plus, 0:cbd:cbd*(sd-1), (1:ad)');
            ind = ind(:)' + (i-1)*ad;
            ind = shift + ind;
            block = self.U(:,ind)'*M*self.U(:,ind);
        end
        
        function block = blockOfGroupAlgebraMatrix(self, M, r, c)
            if nargin < 4 || c == 0
                block = self.blockOfGroupAlgebraMatrix(M, r, 1);
                for c = 2:self.multiplicities(r)
                    block = block + self.blockOfGroupAlgebraMatrix(M, r, c);
                end
                block = block / self.multiplicities(r);
                return
            end
            sizes = self.multiplicities.*self.dimensions;
            shift = sum(sizes(1:r-1)) + (c-1)*self.dimensions(r);
            ind = shift+(1:self.dimensions(r));
            block = self.U(:, ind)'*M*self.U(:, ind);
        end

        function blocks = blocksOfGroupAlgebraMatrix(self, M, applyReynolds)
        % Returns the blocks of a matrix M that represents the group algebra
        %
        % applyReynolds: whether to average M over copies
        %                (optional, default = true)
            if nargin < 3 || applyReynolds
                c = 0;
            else
                c = 1;
            end
            blocks = cell(1, self.nComponents);
            for r = 1:self.nComponents
                blocks{r} = self.blockOfGroupAlgebraMatrix(M, r, c);
            end
        end

        function blocks = blocksOfCentralizer(self, M, applyReynolds)
        % Returns the blocks of a matrix M that commutes with the representation
        % after block diagonalization.
        %
        % applyReynolds: whether to project M in the centralizer algebra first
        %                (optional, default = true)
            if nargin < 3 || applyReynolds
                c = 0;
            else
                c = 1;
            end
            blocks = cell(1, self.nComponents);
            for r = 1:self.nComponents
                blocks{r} = self.blockOfCentralizer(M, r, c);
            end
        end
        
    end
    
    methods (Static)
        
        function irr = fromIsotypicDecomposition(id)
            irrDec = replab.rep.IrrDec.fromIsoDec(id.isoDec);
            U = irrDec.U;
            n = irrDec.nComponents;
            divisionAlgebras = cell(1, n);
            for i = 1:n
                switch irrDec.repTypes(i)
                  case 1 % real
                      divisionAlgebras{i} = replab.rep.DivisionAlgebra.real;
                  case 2 % complex
                      divisionAlgebras{i} = replab.rep.DivisionAlgebra.complex;
                  case 3 % quaternion
                      divisionAlgebras{i} = replab.rep.DivisionAlgebra.quaternion;
                end
            end
            irr = replab.IrreducibleDecomposition(id.rep, U, irrDec.repDims, irrDec.repMuls, divisionAlgebras);
        end
        
    end
end
