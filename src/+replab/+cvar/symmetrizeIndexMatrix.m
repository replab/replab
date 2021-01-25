function indexMatrix = symmetrizeIndexMatrix(indexMatrix, generators, matrixType, field)
% Imposes invariance on a matrices of indices
%
% An index matrix is a matrix with integer entries. The value of every
% element refers to the name of the variable that is meant to appear in
% this entry. An index matrix thus describes a matrix structure.
%
% When requesting that a matrix with given structure be invariant under
% a joint permutation of rows and columns, several entries can typically be
% identifies. The number of variables can then be reduces, as well as the 
%
% Indices values that appear in the matrix are not meaningful, only the
% fact that two indices are identical or different matters.
%
% For convenience, the input is allowed to have negative values as well.
%
% Args:
%   indexMatrix (integer (\*,\*)): matrix of indices
%   generators (permutation): a list of generators under which
%     the matrix remains unchanged
%   matrixType (charstring): one of the following:
%     'full' : no particular structure
%     'symmetric' : transpose-invariant matrix
%     'antisymmetric' : antisymmetric matrix
%
% Returns:
%   integer (\*,\*): symmetrizes matrix of indices
%
% Example:
%   >>> indexMatrix = [1 1 3 4; 1 5 6 30; 3 6 10 11; 4 30 11 15];
%   >>> symIndexMatrix = replab.cvar.symmetrizeIndexMatrix(indexMatrix, {[4 1 2 3]}, 'symmetric');

    % Basic tests
    assert(max(max(abs(indexMatrix - round(indexMatrix)))) == 0, 'The indexMatrix must be a matrix of integers.');
    d = size(indexMatrix, 1);
    for i = 1:length(generators)
        assert(length(generators{i}) == d, 'Generators and indexMatrix dimensions don''t match.');
    end
    assert(isequal(matrixType, 'full') || isequal(matrixType, 'symmetric'), 'The matrix type must be ''full'', ''symmetric''.');
    assert(min(min(indexMatrix)) >= 1, 'Only positive indices are currently supported');
    
    % We renumber all indices so they match default numbering
    [values, indices] = unique(indexMatrix(:), 'first');
    invPerm = sparse(values, 1, indices);
    indexMatrix = full(invPerm(indexMatrix));

    % We write the action of the generators on the matrix elements
    generators2 = cell(size(generators));
    M = reshape(1:d^2, [d d]);
    for i = 1:length(generators)
        generators2{i} = reshape(M(generators{i}, generators{i}), 1, d^2);
    end
    group = replab.PermutationGroup.of(generators2{:});

    % Identify the orbits
    orbits = group.orbits.blockIndex;

    % Make sure orbits are numbered in a similar way
    [values, indices] = unique(orbits(:), 'first');
    invPerm = sparse(values, 1, indices);
    orbits = full(invPerm(orbits));

    % Merge orbits with indices
    pairs = [reshape(indexMatrix, d^2, 1), orbits];

    if isequal(matrixType, 'symmetric')
        % Also impose that indexMatrix is symmetric
        pairs = [pairs; reshape(indexMatrix, d^2, 1) reshape(indexMatrix', d^2, 1)];
    end
    pairs = unique(sort(pairs, 2), 'rows');

    % We identify all connected subsets and attribute a distinct number to
    % their respective elements
    images = replab.graph.connectedComponents(max(pairs, [], 'all'), pairs).';

    % First we identify the index for each element
    % We substitute just one index per class of indices
    [values, indices] = unique(unique(indexMatrix(:)), 'first');
    invPerm = sparse(values, 1, indices);
    tmp = sparse(1:numel(indexMatrix), invPerm(reshape(indexMatrix,1,numel(indexMatrix))), true);
    
    indexMatrix = reshape(tmp*images(values), d, d);
end
