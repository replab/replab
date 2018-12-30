function rid = irreducibleDecomposition(realRep)
% Computes and return the irreducible decomposition of a real representation
    iso = replab.rep.IsoDec.fromAlgebra(realRep.centralizerAlgebra);
    I = replab.rep.IrrDec.fromIsoDec(iso);
    nC = I.nComponents;
    nG = realRep.group.nGenerators;
    components = cell(1, nC);
    shift = 0;
    foundTrivial = false;
    for c = 1:nC
        d = I.repDims(c);
        m = I.repMuls(c);
        switch I.repTypes(c)
          case 1
            da = replab.DivisionAlgebra.real;
          case 2
            da = replab.DivisionAlgebra.complex;
          case 3
            da = replab.DivisionAlgebra.quaternion;
          otherwise
            error('Unknown type');
        end
        ind = shift + (1:m*d);        
        components{c} = replab.RealIrrep.fromParentRealRep(realRep, d, m, da, I.U(:,ind), I.U(:,ind)');
        if components{c}.isTrivial
            % Put the trivial representation first
            assert(~foundTrivial);
            tmp = components{1};
            components{1} = components{c};
            components{c} = tmp;
            foundTrivial = true;
        end
        shift = shift + m*d;
    end
    rid = replab.RealDecompositionRep(realRep, components);
end
