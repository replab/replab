classdef IrreducibleDecomposition < replab.Str
% Irreducible decomposition of a representation

    properties (SetAccess = protected)
        rep; % Representation being decomposed
        centralizerAlgebraBlocks; % Block diagonal matrix algebra
        groupAlgebraBlocks; % Block diagonal matrix algebra
        U;
        n; % Number of isotypic components
        d; % 1xn vector: dimension of representations
        m; % 1xn vector: multiplicity of representations
        divisionAlgebras;
    end

    methods
        
        function self = IrreducibleDecomposition(rep, U, d, m, divisionAlgebras)
            self.rep = rep;
            blockFun = @(n) replab.Matrices(n, n, false, false);
            self.centralizerAlgebraBlocks = arrayfun(blockFun, m, 'UniformOutput', false);
            self.groupAlgebraBlocks = arrayfun(blockFun, d, 'UniformOutput', false);
            self.U = U;
            self.n = length(d);
            self.d = d;
            self.m = m;
            self.divisionAlgebras = divisionAlgebras;
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
