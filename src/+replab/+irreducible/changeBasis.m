function newIrrep2 = changeBasis(irrep1, irrep2, F, G)
% Given two equivalent irreducible subrepresentations, updates the basis of the second one to make them identical
%
% Args:
%   irrep1 (`+replab.SubRep`): First irreducible subrepresentation
%   irrep2 (`+replab.SubRep`): Second irreducible subrepresentation
%   F (double(\*,\*)): Equivariant sample ``irrep1.image(g) * F == F * irrep2.image(g)``
%   G (double(\*,\*)): Equivariant sample ``irrep2.image(g) * G == G * irrep1.image(g)``
    g = irrep1.group.sample;
    tol = replab.globals.doubleEigTol;
    assert(norm(irrep1.image(g) * F - F * irrep2.image(g), 'fro') <= tol);
    assert(norm(irrep2.image(g) * G - G * irrep1.image(g), 'fro') <= tol);
    d = irrep1.dimension;
    assert(irrep2.dimension == d);
    t = trace(F*G);
    F = F * sqrt(abs(d/t));
    G = sign(t) * G * sqrt(abs(d/t));
    FG = F*G;
    norm(FG-eye(d),'fro')
    S = eye(d); % compute the inverse matrix square root of FG
    for i = 1:6
        S = 2*S/(eye(d)+FG*S*S);
    end
    F = S*F;
    G = G*S;
    norm(F*G-eye(d),'fro')
    newI = irrep2.injection*G;
    newP = F*irrep2.projection;
    newIrrep2 = irrep2.parent.subRep(newI, 'projection', newP, 'isIrreducible', true, 'frobeniusSchurIndicator', irrep2.frobeniusSchurIndicator, 'divisionAlgebraName', irrep2.divisionAlgebraName, 'isUnitary', irrep2.isUnitary);
end
