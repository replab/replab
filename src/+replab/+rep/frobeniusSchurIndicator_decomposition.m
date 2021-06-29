function f = frobeniusSchurIndicator_decomposition(rep)
% Computes the Frobenius-Schur indicator by decomposing a representation into irreps
%
% Args:
%   rep (`+replab.Rep`): Representation to compute the Frobenius-Schur indicator of
%
% Returns:
%   integer: Frobenius-Schur indicator
    dec = rep.decomposition;
    f = 0;
    for i = 1:dec.nComponents
        c = dec.component(i);
        if c.nIrreps > 0
            f = f + c.multiplicity * c.irrep(1).frobeniusSchurIndicator;
        end
    end
end
