classdef IrreducibleCommutant < replab.Commutant

    properties (SetAccess = protected)
        irreducible % replab.Irreducible: Decomposition of the representation into irreducibles
    end
    
    methods (Static, Access = protected)
        
        function block = extractSimpleBlockInBasis(X, U, shift, cd, m)
        % Changes the basis and projects a block from a commutant element (real-type real, or complex)
        %
        % We have that self.irreducible.parent.image(g) * X == X * self.irreducible.parent.image(g) for all g in the group.
        %
        % Args:
        %   X (double): Matrix that commutes with the irreducible representation, in the original representation basis
        %   U (double): Change of basis such that U * X * U' is block diagonal in the irreducible decomposition basis
        %   shift (integer): Cumulated dimension of all previous isotypic components
        %   cd (integer): Dimension of a single copy of the irreducible representation in the isotypic component
        %   m (integer): Multiplicity of the irreducible representation in the isotypic component
        %
        % Returns:
        %   double: Block corresponding to the isotypic component
            block = zeros(m, m);
            for i = 1:cd
                Ublock = U(shift+(i:cd:m*cd), :);
                block = block + Ublock * X * Ublock';
            end
            block = block/cd;
        end
        
        function block = extractSimpleBlock(X, shift, cd, m)
        % Projects a block from a commutant element (real-type real, or complex)
        %
        % We have that self.irreducibleRep.image(g) * X == X * self.irreducibleRep.image(g) for all g in the group.
        %
        % Args:
        %   X (double): Matrix that commutes with the irreducible representation, in the block-diagonalization basis
        %   shift (integer): Cumulated dimension of all previous isotypic components
        %   cd (integer): Dimension of a single copy of the irreducible representation in the isotypic component
        %   m (integer): Multiplicity of the irreducible representation in the isotypic component
        %
        % Returns:
        %   double: Block corresponding to the isotypic component
            block = zeros(m, m);
            for i = 1:cd
                block = block + X(shift+(i:cd:m*cd), shift+(i:cd:m*cd));
            end
            block = block/cd;
        end
        
        function [A B] = extractComplexTypeBlockInBasis(X, U, shift, cd, m)
        % Changes the basis and projects a block from a commutant element (complex-type real)
        %
        % See replab.IrreducibleCommutant.extractRealBlockInBasis for arguments
            A = zeros(m, m);
            B = zeros(m, m);
            for i = 1:2:cd
                U1 = U(shift+(i:cd:m*cd), :);
                U2 = U(shift+(i:cd:m*cd)+1, :);
                A = A + U1*X*U1' + U2*X*U2';
                B = B + U2*X*U1' - U1*X*U2';
            end
            A = A/cd;
            B = B/cd;
        end
        
        function [A B] = extractComplexTypeBlock(X, shift, cd, m)
        % Projects a block from a commutant element (complex-type real)
        %
        % See replab.IrreducibleCommutant.extractRealBlock for arguments
            A = zeros(m, m);
            B = zeros(m, m);
            for i = 1:2:cd
                r = shift+(i:cd:m*cd);
                A = A + X(r, r) + X(r+1, r+1);
                B = B + X(r+1, r) - X(r, r+1);
            end
            A = A/cd;
            B = B/cd;
        end
        
        function [A B C D] = extractQuaternionTypeBlock(X, shift, cd, m)
        % Projects a block from a commutant element (quaternion-type real)
        %
        % See replab.IrreducibleCommutant.extractRealBlock for arguments
            A = zeros(m, m);
            B = zeros(m, m);
            C = zeros(m, m);
            D = zeros(m, m);
            % shape of things that commute with our representation quaternion encoding
            % [ a -b -c -d
            %   b  a  d -c
            %   c -d  a  b
            %   d  c -b  a]
            for i = 1:4:cd
                r1 = shift+(i:cd:m*cd);
                r2 = shift+(i:cd:m*cd)+1;
                r3 = shift+(i:cd:m*cd)+2;
                r4 = shift+(i:cd:m*cd)+3;
                A = A + X(r1,r1) + X(r2,r2) + X(r3,r3) + X(r4,r4);
                B = B + X(r2,r1) - X(r1,r2) - X(r4,r3) + X(r3,r4);
                C = C + X(r3,r1) - X(r1,r3) - X(r2,r4) + X(r4,r2);
                D = D + X(r4,r1) - X(r3,r2) + X(r2,r3) - X(r1,r4);
            end
            A = A/cd;
            B = B/cd;
            C = C/cd;
            D = D/cd;
        end
        
        function [A B C D] = extractQuaternionTypeBlockInBasis(X, U, shift, cd, m)
        % Changes the basis and projects a block from a commutant element (quaternion-type real)
        %
        % See replab.IrreducibleCommutant.extractRealBlockInBasis for arguments
            A = zeros(m, m);
            B = zeros(m, m);
            C = zeros(m, m);
            D = zeros(m, m);
            % shape of things that commute with our representation quaternion encoding
            % [ a -b -c -d
            %   b  a  d -c
            %   c -d  a  b
            %   d  c -b  a]
            for i = 1:4:cd
                U1 = U(shift+(i:cd:m*cd), :);
                U2 = U(shift+(i:cd:m*cd)+1, :);
                U3 = U(shift+(i:cd:m*cd)+2, :);
                U4 = U(shift+(i:cd:m*cd)+3, :);
                A = A + U1*X*U1' + U2*X*U2' + U3*X*U3' + U4*X*U4';
                B = B + U2*X*U1' - U1*X*U2' - U4*X*U3' + U3*X*U4';
                C = C + U3*X*U1' - U1*X*U3' - U2*X*U4' + U4*X*U2';
                D = D + U4*X*U1' - U3*X*U2' + U2*X*U3' - U1*X*U4';
            end
            A = A/cd;
            B = B/cd;
            C = C/cd;
            D = D/cd;
        end
        
    
    end
    
    methods
        
        function self = IrreducibleCommutant(irreducibleRep)
        % Constructs the commutant of a irreducible decomposition representation
        %
        % Args:
        %   irreducibleRep (replab.IrreducibleRep): Irreducible decomposition representation
            assert(isa(irreducibleRep, 'replab.IrreducibleRep'));
            self = self@replab.Commutant(irreducibleRep);
            self.irreducible = irreducibleRep.irreducible;
        end
        
        function s = singleBlockSize(self, i)
        % Returns the size of a commutant algebra element block, without repetition
            iso = self.irreducible.component(i);
            d = iso.dimension;
            m = iso.multiplicity;
            cd = iso.copyDimension;
            first = iso.copy(1);
            if first.overC || isequal(first.irrepInfo.divisionAlgebra, 'R')
                s = m;
            else
                switch first.irrepInfo.divisionAlgebra
                  case 'C'
                    s = m*2;
                  case 'H'
                    s = m*4;
                end
            end
        end
        
        function block = blockOfParentElement(self, X, i, withRepetitions)
        % Extracts a block from an element of the commutant algebra expressed in the original basis
        %
        % Args:
        %   X (double): Matrix that commutes with the original representation `self.irreducible.rep`
        %   i (integer): Index of the block to extract
        %   withRepetitions (logical): Whether to reproduce the degeneracy due to the representation dimension
        %
        % Returns:
        %   double: The `i`-th block after change of basis and projection
            dimensions = cellfun(@(x) x.dimension, self.irreducible.components);
            shift = sum(dimensions(1:i-1));
            iso = self.irreducible.component(i);
            first = iso.copy(1);
            d = iso.dimension;
            m = iso.multiplicity;
            cd = iso.copyDimension;
            U = self.irreducible.U;
            if first.overC || isequal(first.irrepInfo.divisionAlgebra, 'R')
                % Complex representation or real-type real representation
                block = replab.IrreducibleCommutant.extractSimpleBlockInBasis(X, U, shift, cd, m);
                if withRepetitions
                    block = kron(block, eye(cd));
                end
            else
                switch first.irrepInfo.divisionAlgebra
                  case 'C'
                    [A B] = replab.IrreducibleCommutant.extractComplexTypeBlockInBasis(X, U, shift, cd, m);
                    if withRepetitions
                        block = kron(A, eye(cd)) + kron(B, kron(eye(cd/2), [0 -1; 1 0]));
                    else
                        block = kron(A, eye(2)) + kron(B, [0 -1; 1 0]);
                    end
                  case 'H'
                    % Quaternion-type real representation
                    [A B C D] = replab.IrreducibleCommutant.extractQuaternionTypeBlockInBasis(X, U, shift, cd, m);
                    basisB = [ 0 -1  0  0
                               1  0  0  0
                               0  0  0  1
                               0  0 -1  0];
                    basisC = [ 0  0 -1  0
                               0  0  0 -1
                               1  0  0  0
                               0  1  0  0];
                    basisD = [ 0  0  0 -1
                               0  0  1  0
                               0 -1  0  0
                               1  0  0  0];
                    if withRepetitions
                    block = kron(A, eye(cd)) + kron(B, kron(eye(cd/4), basisB)) + ...
                            kron(C, kron(eye(cd/4), basisC)) + kron(D, kron(eye(cd/4), basisD));
                    else
                        block = kron(A, eye(4)) + kron(B, basisB) + kron(C, basisC) + kron(D, basisD);
                    end
                end
            end
        end
        
        function block = projectedBlock(self, X, blockIndex, withRepetitions)
        % Returns a block of a matrix projected in the commutant algebra
        %
        % Args:
        %   X (double): Matrix to project on the commutant algebra
        %   blockIndex (integer): Isotypic component index
        %   withRepetitions (logical): Whether to reproduce the degeneracy due to the representation dimension
        %
        % Returns:
        %   double: The projected block
            dimensions = cellfun(@(x) x.dimension, self.irreducible.components);
            shift = sum(dimensions(1:blockIndex-1));
            iso = self.irreducible.component(blockIndex);
            first = iso.copy(1);
            d = iso.dimension;
            m = iso.multiplicity;
            cd = iso.copyDimension;
            if first.overC || isequal(first.irrepInfo.divisionAlgebra, 'R')
                % Complex representation or real-type real representation
                block = replab.IrreducibleCommutant.extractSimpleBlock(X, shift, cd, m);
                if withRepetitions
                    block = kron(block, eye(cd));
                end
            else
                switch first.irrepInfo.divisionAlgebra
                  case 'C'
                    [A B] = replab.IrreducibleCommutant.extractComplexTypeBlock(X, shift, cd, m);
                    if withRepetitions
                        block = kron(A, eye(cd)) + kron(B, kron(eye(cd/2), [0 -1; 1 0]));
                    else
                        block = kron(A, eye(2)) + kron(B, [0 -1; 1 0]);
                    end
                  case 'H'
                    % Quaternion-type real representation
                    [A B C D] = replab.IrreducibleCommutant.extractQuaternionTypeBlock(X, shift, cd, m);
                    basisB = [ 0 -1  0  0
                               1  0  0  0
                               0  0  0  1
                               0  0 -1  0];
                    basisC = [ 0  0 -1  0
                               0  0  0 -1
                               1  0  0  0
                               0  1  0  0];
                    basisD = [ 0  0  0 -1
                               0  0  1  0
                               0 -1  0  0
                               1  0  0  0];
                    if withRepetitions
                    block = kron(A, eye(cd)) + kron(B, kron(eye(cd/4), basisB)) + ...
                            kron(C, kron(eye(cd/4), basisC)) + kron(D, kron(eye(cd/4), basisD));
                    else
                        block = kron(A, eye(4)) + kron(B, basisB) + kron(C, basisC) + kron(D, basisD);
                    end
                end
            end
        end
        
        function X1 = project(self, X)
            I = self.irreducible;
            blocks = cell(1, I.nComponents);
            for c = 1:I.nComponents
                blocks{c} = self.projectedBlock(X, c, true);
            end
            X1 = blkdiag(blocks{:});
        end
        
    end
    
end
