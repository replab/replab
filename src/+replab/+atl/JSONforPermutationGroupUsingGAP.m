function json = JSONforPermutationGroupUsingGAP(G, irreps)
% Runs GAP System to compute the character table/representation information about a permutation group
%
% Args:
%   G (`+replab.PermutationGroup`): Group to compute the information of
%   irreps (logical): Whether to compute irreducible representations
%
% Returns:
%   charstring: The JSON data corresponding to the group
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, replab.atl.getGAPScript(G, irreps));
    fclose(fid);
    [status, json] = system([replab.globals.gapBinaryPath ' -n -q <' tfile]);
    ind = find(json == '{', 1);
    json = json(ind:end);
    delete(tfile);
end
