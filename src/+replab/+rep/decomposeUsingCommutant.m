function sub = decomposeUsingCommutant(rep)
% Decomposes the given representation into irreducible representations
%
% Uses eigenvalue decomposition on a generic commutant sample
%
% Args:
%   rep (replab.Rep): The representation to decompose
%
% Returns:
%   row cell array of replab.SubRep: List of subrepresentations
    assert(isa(rep, 'replab.Rep'));
    tol = replab.Settings.doubleEigTol;
    % Sample a self adjoint commutant element 
    C = full(rep.commutant.sampleSelfAdjoint);
    [U D] = replab.rep.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    sub = cell(1, n);
    for i = 1:n
        sub{i} = rep.subRepUnitary(U(:, runs{i})');
    end
end
