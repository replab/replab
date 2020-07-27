function basis = findIsotypicBasis(irrep, rep)
% Returns a basis of an isotypic component, with the multiplicity space identified
%
% Args:
%   irrep (`+replab.Rep`): Irreducible representatoin
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   double(\*,\*), may be sparse: Basis of the isotypic component corresponding to ``irrep`` of dimension ``rep.dimension x (multiplicity*irrep.dimension)``
end
