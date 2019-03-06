% Here we explore some matrix structures

%generators = {[2 3 4 6 5 8 7 1]}
generators = {[2 3 1]}
%generators = {[2 3 1 5 6 4]}

n = size(generators{1},2)

group = replab.Permutations(n).subgroup(generators)

irrDecomp = group.naturalRepresentation.irreducible

U = irrDecomp.U


dimensions1 = zeros(1, irrDecomp.nComponents);
multiplicities = zeros(1, irrDecomp.nComponents);
types = '';
for i = 1:irrDecomp.nComponents
    dimensions1(i) = irrDecomp.component(i).dimension1;
    multiplicities(i) = irrDecomp.component(i).multiplicity;
    types(i) = irrDecomp.component(i).divisionAlgebra.shortName;
end
dimensions1
multiplicities
types


rho = group.naturalRepresentation

M = group.naturalRepresentation.centralizerAlgebra.project(rand(n,n))

U'*M*U

