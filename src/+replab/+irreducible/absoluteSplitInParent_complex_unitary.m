function subs = absoluteSplitInParent_complex_unitary(sub, sample)
% Decomposes a complex unitary subrepresentation into subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Unitary subrepresentation to split further
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%x
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.isUnitary);
    assert(sub.overC);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    subI = sub.injection('double/sparse');
    subP = sub.projection('double/sparse');
    S = full(subP * sample * subI);
    X = (S + S')/2;
% $$$     if isa(sub.parent.commutant, 'replab.equi.Equivariant_forCompactGroup') && replab.globals.useReconstruction
% $$$         B = sub.parent.commutant.blocks;
% $$$         s = [];
% $$$         for i = 1:size(B, 2)
% $$$             if i ~= 1
% $$$                 s = [s '+'];
% $$$             end
% $$$             s = [s sprintf('%d', length(B{1,i}))];
% $$$         end
% $$$         replab.msg(1, 'Blocks of size %s', s);
% $$$         U = zeros(sub.parent.dimension, sub.parent.dimension);
% $$$         D = U;
% $$$         shift = 0;
% $$$         for i = 1:size(B, 2)
% $$$             b1 = B{1,i};
% $$$             b2 = B{2,i};
% $$$             assert(isequal(sort(b1), sort(b2)));
% $$$             [Ub, Db] = eig(X(b1, b1));
% $$$             l = length(b1);
% $$$             U(b1, shift+(1:l)) = Ub;
% $$$             D(shift+(1:l), shift+(1:l)) = Db;
% $$$             shift = shift + l;
% $$$         end
% $$$     else
    [U, D] = eig(X);
        %end
    D = reshape(diag(D), [1 d]);
    P = replab.Partition.fromApproximateVector(D, tol);
    blocks = P.blocks;
    n = P.nBlocks;
    if n == 1
        sub.cache('isIrreducible', true, '==');
        subs = {sub};
    else
        subs = cell(1, n);
        for i = 1:n
            basis = U(:, blocks{i});
            I = subI * basis;
            if sub.mapsAreAdjoint
                P = I';
            else
                P = basis' * subP;
            end
            subs{i} = sub.parent.subRep(I, 'projection', P, 'isUnitary', true);
        end
    end
end
