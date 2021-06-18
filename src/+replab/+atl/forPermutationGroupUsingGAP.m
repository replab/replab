function result = forPermutationGroupUsingGAP(G, irreps)
% Runs GAP System to compute the character table/representation information about a permutation group
%
% Args:
%   G (`.PermutationGroup`): Group to compute the information of
%   irreps (logical): Whether to compute irreducible representations
%
% Returns:
%   `+replab.AbstractGroup`: The computed abstract group with character table
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, replab.atl.getGAPScript(G, irreps));
    fclose(fid);
    [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
    delete(tfile);
    %A = replab.AtlasEntry.groupFromJSONData(replab.util.parseJSON(result));
end
