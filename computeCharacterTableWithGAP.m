function computeCharacterTableWithGAP(group, computeIrreps)
% Populates the character table of a finite group
%
% Args:
%   group (`.FiniteGroup`): Group without complex character table
%   computeIrreps (logical, optional): Whether to compute complex irreducible representations, default: true
    if nargin < 2 || isempty(computeIrreps)
        computeIrreps = true;
    end
    if isa(group, 'replab.PermutationGroup')
        permGroup = group;
    else
        permGroup = group.niceGroup;
    end
    tfile = tempname();
    fid = fopen(tfile, 'wt');
    fprintf(fid, 'G := Group(%s);;', strjoin(cellfun(@(g) replab.atl.permToGap(g), group.generators, 'uniform', 0), ', '));
    fwrite(fid, fileread(fullfile(replab.globals.replabPath, 'computeCharacterTableWithGAP.g')));
    fclose(fid);
    [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
    delete(tfile);
    data = replab.util.parseJSON(result);
    irreps = replab.atl.parseIrreps(group, 'C', group.generators, data);
    C = replab.ComplexCharacterTable.fromIrreps(group, irreps);
    group.cache('complexCharacterTable', C, 'error');
end
