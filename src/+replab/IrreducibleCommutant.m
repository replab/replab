classdef IrreducibleCommutant < replab.Commutant

    properties (SetAccess = protected)
        irreducible % replab.Irreducible: Decomposition of the representation into irreducibles
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
        
        function X1 = project(self, X)
            I = self.irreducible;
            blocks = {};
            shift = 0;
            for c = 1:I.nComponents
                iso = self.irreducible.component(c);
                first = iso.copy(1);
                d = iso.dimension;
                m = iso.multiplicity;
                cd = iso.copyDimension;
                if self.irreducible.parent.overC || isequal(first.irrepInfo.divisionAlgebra, 'R')
                    block = zeros(m, m);
                    for i = 1:cd
                        block = block + X(shift+(i:cd:d), shift+(i:cd:d));
                    end
                    block = block/cd;
                    block = kron(block, eye(cd));
                elseif isequal(first.irrepInfo.divisionAlgebra, 'C')
                    A = zeros(m, m);
                    B = zeros(m, m);
                    for i = 1:2:cd
                        A = A + X(shift+(i:cd:d), shift+(i:cd:d)) + X(shift+(i+1:cd:d), shift+(i+1:cd:d));
                        B = B + X(shift+(i+1:cd:d), shift+(i:cd:d)) - X(shift+(i:cd:d), shift+(i+1:cd:d));
                    end
                    A = A/cd;
                    B = B/cd;
                    block = kron(A, eye(cd)) + kron(B, kron(eye(cd/2), [0 -1; 1 0]));
                else
                    assert(isequal(first.irrepInfo.divisionAlgebra, 'H'));
                    A = zeros(m, m);
                    B = zeros(m, m);
                    C = zeros(m, m);
                    D = zeros(m, m);
                    for i = 1:4:cd
                        A = A + X(shift+(i:cd:d),shift+(i:cd:d)) + X(shift+(i+1:cd:d),shift+(i+1:cd:d)) ...
                            + X(shift+(i+2:cd:d),shift+(i+2:cd:d)) + X(shift+(i+3:cd:d),shift+(i+3:cd:d));
                        B = B + X(shift+(i+1:cd:d),shift+(i:cd:d)) - X(shift+(i:cd:d),shift+(i+1:cd:d)) ...
                            + X(shift+(i+3:cd:d),shift+(i+2:cd:d)) - X(shift+(i+2:cd:d),shift+(i+3:cd:d));
                        C = C + X(shift+(i+2:cd:d),shift+(i:cd:d)) - X(shift+(i:cd:d),shift+(i+2:cd:d)) ...
                            + X(shift+(i+1:cd:d),shift+(i+3:cd:d)) - X(shift+(i+3:cd:d),shift+(i+1:cd:d));
                        D = D + X(shift+(i+3:cd:d),shift+(i:cd:d)) + X(shift+(i+2:cd:d),shift+(i+1:cd:d)) ...
                            - X(shift+(i+1:cd:d),shift+(i+2:cd:d)) - X(shift+(i:cd:d),shift+(i+3:cd:d));
                    end
                    A = A/cd;
                    B = B/cd;
                    C = C/cd;
                    D = D/cd;
                    basisB = [ 0 -1  0  0
                               1  0  0  0
                               0  0  0 -1
                               0  0  1  0];
                    basisC = [ 0  0 -1  0
                               0  0  0  1
                               1  0  0  0
                               0 -1  0  0];
                    basisD = [ 0  0  0 -1
                               0  0 -1  0
                               0  1  0  0
                               1  0  0  0];
                    block = kron(A, eye(cd)) + kron(B, kron(eye(cd/4), basisB)) + ...
                            kron(C, kron(eye(cd/4), basisC)) + kron(D, kron(eye(cd/4), basisD));
                end
                blocks{c} = block;
                shift = shift + d;
            end
            X1 = blkdiag(blocks{:});
        end
        
    end
    
end
