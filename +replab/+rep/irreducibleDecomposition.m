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
        images1 = cell(1, nG);
        imagesInv1 = cell(1, nG);
        for i = 1:nG
            rho = realRep.images{i};
            rhoInv = realRep.imagesInv{i};
            img = zeros(d, d);
            imgInv = zeros(d, d);
            for j = 1:m
                ind = shift + (j-1)*d + (1:d);
                img = img + I.U(:,ind)'*rho*I.U(:,ind);
                imgInv = imgInv + I.U(:,ind)'*rhoInv*I.U(:,ind);
            end
            img = img / m;
            imgInv = imgInv / m;
            img = da.projectMatrix(img);
            imgInv = da.projectMatrix(imgInv);
            images1{i} = img;
            imagesInv1{i} = imgInv;
        end
        ind = shift + (1:m*d);
        components{c} = replab.RealIrrep(realRep.group, d, m, da, images1, imagesInv1, ...
                                         realRep, I.U(:,ind), I.U(:,ind)');
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
