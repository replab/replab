function A = forPermutationGroupUsingGAP(G, irreps)
% Runs GAP System to compute the character table/representation information about a permutation group
%
% Args:
%   G (`.PermutationGroup`): Group to compute the information of
%   irreps (logical): Whether to compute irreducible representations
%
% Returns:
%   `.AtlasEntry`: The completed atlas for the given group
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, replab.AtlasEntry.getGAPScript(G, irreps));
    fclose(fid);
    [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
    delete(tfile);
    A = replab.AtlasEntry.parse(result);
end
