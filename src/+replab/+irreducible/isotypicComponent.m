function iso = isotypicComponent(rep, irrep, type)
% Returns the isotypic component corresponding to a given irrep
%
% Based on original code by Sotiris Mygdalas, summer 2020
%
% Args:
%   rep (`+replab.Rep`): Representation to find the isotypic component of
%   irrep (`+replab.Rep`): Irreducible representation
%   type ('exact', 'double', 'double/sparse'): Type of the projector to compute
%
% Returns:
%   `+replab.Isotypic`: Isotypic component, may be empty
    if nargin < 3 || isempty(type)
        type = 'double';
    end
    assert(rep.field == irrep.field);
    if irrep.overR
        assert(irrep.frobeniusSchurIndicator == 1, 'Over the reals, isotypicComponent only support real-type/absolutely irreducible representations');
    end
    D = rep.dimension;
    d = irrep.dimension;
    T = replab.irreducible.serreTensor(rep, irrep, type);
    proj1 = T(:,:,1,1);
    % find a linear basis of the range of proj1
    % we use the row reduced echelon form, which finds indices of linearly independent columns
    [~, pivot] = rref(double(proj1));
    m = length(pivot);
    if strcmp(type, 'exact')
        I = replab.cyclotomic.zeros(D, m, d);
    else
        I = zeros(D, m, d);
    end
    % we create the first element of the injection map
    I(:,:,1) = T(:,pivot,1,1);
    for i = 2:d
        % and compute the others using Serre's trick
        I(:,:,i) = T(:,:,i,1) * T(:,pivot,1,1);
    end
    I = permute(I, [1 3 2]);
    sub = rep.subRep(reshape(I, [D m*d]));
    Piso = sub.projection(type);
    P = permute(reshape(Piso, [d m D]), [1 3 2]);
    irreps = cell(1, m);
    for i = 1:m
        irreps{i} = rep.subRep(I(:,:,i), 'projection', P(:,:,i), 'isIrreducible', true, 'isUnitary', irrep.isUnitary);
    end
    iso = replab.Isotypic(rep, irreps, irrep, Piso);
end
