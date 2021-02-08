function [iso, zeroErrors, nonZeroErrors] = findIsotypic(parent, irreps, sample)
% Identifies equivalent irreducible subrepresentations and sorts them in harmonized isotypic components
%
% All subrepresentations must have the same parent.
%
% Args:
%   parent (`+replab.Rep`): Parent representation of which the irreps have been computed
%   irreps (cell(1,\*) of `+replab.SubRep`): Biorthogonal irreducible subrepresentations with their division algebra in canonical form
%   sample (double(\*,\*)): Sample of ``parent.commutant``
    dims = cellfun(@(r) r.dimension, irreps);
    zeroErrors = zeros(1, 0);
    nonZeroErrors = zeros(1, 0);
    iso = cell(1, 0);
    for d = unique(dims)
        inds = find(dims == d);
        n = length(inds);
        adj = zeros(n, n);
        for i = 1:n
            for j = 1:n
                ri = irreps{inds(i)};
                rj = irreps{inds(j)};
                Eij = ri.projection*sample*rj.injection;
                Eji = rj.projection*sample*ri.injection;
                S = Eij*Eji;
                tol = 1e-10;
                if norm(S) > tol
                    assert(strcmp(ri.divisionAlgebraName, rj.divisionAlgebraName));
                    if strcmp(ri.divisionAlgebraName, 'quaternion.rep')
                        S1 = replab.irreducible.projectScalar(S, 'quaternion.equivariant');
                    else
                        S1 = replab.irreducible.projectScalar(S, ri.divisionAlgebraName);
                    end
                    err = norm(S - S1, 'fro');
                    nonZeroErrors(1,end+1) = err;
                    adj(i,j) = 1;
                    adj(j,i) = 1;
                else
                    err = norm(S, 'fro');
                    zeroErrors(1,end+1) = err;
                end
            end
        end
        G = replab.UndirectedGraph.fromAdjacencyMatrix(adj);
        P = G.connectedComponents;
        for i = 1:P.nBlocks
            blk = P.block(i);
            iso{1,end+1} = replab.Isotypic.fromBiorthogonalIrreps(parent, irreps(inds(blk)), d, false);
        end
    end
end
