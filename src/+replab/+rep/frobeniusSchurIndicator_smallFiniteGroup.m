function f = frobeniusSchurIndicator_smallFiniteGroup(rep)
% Computes the Frobenius-Schur indicator using an explicit sum
%
% Args:
%   rep (`+replab.Rep`): Representation to compute the Frobenius-Schur indicator of
%
% Returns:
%   integer: Frobenius-Schur indicator
    % for a finite group, use conjugacy classes
    C = rep.group.conjugacyClasses.classes;
    n = length(C);
    g2 = cellfun(@(c) rep.group.composeN(c.representative, 2), C, 'uniform', 0);
    factor = cellfun(@(c) rep.group.order/c.nElements, C, 'uniform', 0);
    if rep.isExact
        f = replab.cyclotomic.zeros(1, 1);
        for i = 1:n
            f = f + trace(rep.image(g2{i}, 'exact'))/replab.cyclotomic(factor{i});
        end
        f = double(f);
        assert(isreal(f) && round(f) == f);
        f = round(f);
    else
        if rep.errorBound >= 1
            error('Error on this representation is too big to compute the Frobenius-Schur indicator');
        end
        f = 0;
        for i = 1:n
            f = f + trace(rep.image(g2{i}, 'double/sparse'))/double(factor{i});
        end
        f = round(f);
    end
end
