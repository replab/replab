function filename = writeAtlasForPermutationGroupUsingGAP(G, irreps)
% Returns the GAP script that outputs the JSON data corresponding to the given group
%
% Args:
%   filename (charstring): Path to JSON file to write
%   G (`+replab.PermutationGroup`): Group to compute the information of
%   name (charstring): Group name
%   irreps (logical): Whether to compute irreducible representations
%
% Returns:
%   charstring: GAP System script that computes the relevant data and outputs it in JSON form
    json = replab.atl.JSONforPermutationGroupUsingGAP(G, irreps);
    res = regexp(json,'SmallGroup\((\d+),(\d+)\)','tokens');
    assert(length(res) == 1);
    order = res{1}{1};
    id = res{1}{2};
    filename = sprintf('%s_%s.json', order, id);
    path = fullfile(replab.globals.replabPath, 'atlas', filename);
    fid = fopen(path, 'wt');
    fwrite(fid, json);
    fclose(fid);
end
