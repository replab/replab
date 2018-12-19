classdef IrreducibleDecomposition < replab.Str
% Irreducible decomposition of a representation

    properties (SetAccess = protected)
        rep; % Representation being decomposed
        centralizerAlgebraBlocks; % Block diagonal matrix algebra
        groupAlgebraBlocks; % Block diagonal matrix algebra
        n; % Number of isotypic components
        m; % 1xn vector: multiplicity of representations
        d; % 1xn vector: dimension of representations
        divisionAlgebra;
    end

end
