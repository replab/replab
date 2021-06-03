function writeJSONforPermutationGroupUsingGAP(G, filename, irreps)
% Returns the GAP script that outputs the JSON data corresponding to the given group
%
% Args:
%   G (`.PermutationGroup`): Group to compute the information of
%   filename (charstring): Path to JSON file to write
%   irreps (logical): Whether to compute irreducible representations
%
% Returns:
%   charstring: GAP System script that computes the relevant data and outputs it in JSON form
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, replab.AtlasEntry.getGAPScript(G, irreps));
    fclose(fid);
    [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
    fid = fopen(filename, 'wt');
    fwrite(fid, result);
    fclose(fid);
end
