function [iso, zeroErrors, nonZeroErrors] = findIsotypic(parent, irreps, sample)
% Identifies equivalent irreducible subrepresentations and sorts them in harmonized isotypic components
%
% All subrepresentations must have the same parent.
%
% Args:
%   parent (`+replab.Rep`): Parent representation of which the irreps have been computed
%   irreps (cell(1,\*) of `+replab.SubRep`): Biorthogonal irreducible subrepresentations with their division algebra in canonical form
%   sample (double(\*,\*)): Sample of ``parent.commutant``
%
% Returns
% -------
%   iso: cell(1,\*) of `+replab.Isotypic`
%     Isotypic components
%   zeroErrors: double(1,\*)
%     Errors in the blocks determined to be zero
%   nonZeroErrors: double(1,\*)
%     Residual error in product scalar maps
    dims = cellfun(@(r) r.dimension, irreps);
    zeroErrors = zeros(1, 0);
    nonZeroErrors = zeros(1, 0);
    iso = cell(1, 0);
    for d = unique(dims)
        todo = reshape(find(dims == d), 1, []);
        while ~isempty(todo)
            i = todo(1);
            todo = todo(2:end);
            ri = irreps{i};
            rip_sample = ri.projection*sample;
            sample_rii = sample*ri.injection;
            comp = {irreps{i}};
            for j = todo
                rj = irreps{j};
                Eij = rip_sample*rj.injection;
                Eji = rj.projection*sample_rii;
                S = Eij*Eji;
                tol = 1e-10;
                if norm(S, 'fro') > tol
                    todo = todo(todo ~= j);
                    assert(strcmp(ri.divisionAlgebraName, rj.divisionAlgebraName));
                    if strcmp(ri.divisionAlgebraName, 'H->R:rep')
                        S1 = replab.irreducible.projectScalar(S, 'H->R:equivariant');
                    else
                        S1 = replab.irreducible.projectScalar(S, ri.divisionAlgebraName);
                    end
                    err = norm(S - S1, 'fro');
                    nonZeroErrors(1,end+1) = err;
                    comp{1,end+1} = replab.irreducible.changeBasis(ri, rj, Eij, Eji);
                else
                    err = norm(S, 'fro');
                    zeroErrors(1,end+1) = err;
                end
            end
            iso{1,end+1} = replab.Isotypic.fromIrreps(parent, comp, comp{1}, 'irrepsAreBiorthogonal', true, 'irrepsAreHarmonized', true);
        end
    end
    dims = cellfun(@(c) c.irrepDimension, iso);
    muls = cellfun(@(c) c.multiplicity, iso);
    [~, ind] = sortrows([dims(:) muls(:)]);
    iso = iso(ind);
end
