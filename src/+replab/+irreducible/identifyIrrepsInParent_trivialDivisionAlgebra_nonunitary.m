function irreps = identifyIrrepsInParent_trivialDivisionAlgebra_nonunitary(sub, sample)
% Identifies the irreducible representation(s) directly present in a real or complex subrepresentation
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation with `+replab.Rep.divisionAlgebraName` set to ``''``
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: A singleton cell array containing the argument ``sub`` with `+replab.Rep.isIrreducible` updated to true, or an empty array
    assert(strcmp(sub.divisionAlgebraName, ''));
    d = sub.dimension;
    I = sub.injection('double/sparse');
    P = sub.projection('double/sparse');
    % the 2-norm is dominated by the fro-norm
    S = P*sample*I;
    f = trace(S)/d;
    delta = norm(S - f*eye(d), 'fro');
    tol = replab.globals.doubleEigTol; % TODO: better tolerance test
    if delta <= tol
        sub.cache('isIrreducible', true, '==');
        if sub.overR
            sub.cache('frobeniusSchurIndicator', 1, '==');
        end
        irreps = {sub};
    else
        irreps = {};
    end
end
