function newIrrep2 = changeBasis(irrep1, irrep2, F, G)
% Given two equivalent irreducible subrepresentations, updates the basis of the second one to make them identical
%
% Args:
%   irrep1 (`+replab.SubRep`): First irreducible subrepresentation
%   irrep2 (`+replab.SubRep`): Second irreducible subrepresentation
%   F (double(\*,\*)): Equivariant sample ``irrep1.image(g) * F == F * irrep2.image(g)``
%   G (double(\*,\*)): Equivariant sample ``irrep2.image(g) * G == G * irrep1.image(g)``
    d = irrep1.dimension;
    assert(irrep2.dimension == d);
    % force F*G ~ eye(n)
    t = trace(F*G);
    F = F * sqrt(abs(d/t));
    G = sign(t) * G * sqrt(abs(d/t));
    % now approxiate F*G = eye(n) better
    FG = F*G;
    % Compute the inverse matrix square root of F*G
    % We use the Newton-type iteration given in:
    % N. Sherif, "On the computation of a matrix inverse square root"
    % Computing, vol. 46, no. 4, pp. 295–305, Dec. 1991, doi: 10.1007/BF02257775.
    nIterations = 6;
    S = eye(d);
    for i = 1:nIterations
        S = 2*S/(eye(d)+FG*S*S);
    end
    F = S*F;
    G = G*S;
    newI = irrep2.injection*G;
    newP = F*irrep2.projection;
    newIrrep2 = irrep2.parent.subRep(newI, 'projection', newP, 'isIrreducible', true, 'frobeniusSchurIndicator', irrep2.frobeniusSchurIndicator, 'divisionAlgebraName', irrep2.divisionAlgebraName, 'isUnitary', irrep2.isUnitary);
end
