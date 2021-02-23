function X1 = projectScalar(X, divisionAlgebraName)
    D = size(X, 1);
    if isempty(divisionAlgebraName)
        X1 = trace(X)/D*speye(D);
    else
        if strcmp(divisionAlgebraName, 'complex')
            d = 2;
        else
            assert(strcmp(divisionAlgebraName, 'quaternion.equivariant'));
            d = 4;
        end
        blk = X(1:d,1:d);
        for i = 2:D/d
            blk = blk + X((i-1)*d+(1:d),(i-1)*d+(1:d));
        end
        blk = d*blk/D;
        da = replab.DivisionAlgebra(divisionAlgebraName);
        X1 = kron(eye(D/d), da.project(blk));
    end
end
