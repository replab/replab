        function S = getGAPScript(G, irreps)
        % Returns the GAP script that outputs the JSON data corresponding to the given group
        %
        % Args:
        %   G (`.PermutationGroup`): Group to compute the information of
        %   irreps (logical): Whether to compute irreducible representations
        %
        % Returns:
        %   charstring: GAP System script that computes the relevant data and outputs it in JSON form
            if irreps
                scriptName = 'AtlasEntry_irreps.g';
            else
                scriptName = 'AtlasEntry_noirreps.g';
            end
            firstLine = sprintf('G := Group(%s);;', strjoin(cellfun(@(g) replab.AtlasEntry.permToGap(g), G.generators, 'uniform', 0), ', '));
            rest = fileread(fullfile(replab.globals.replabPath, 'src', '+replab', scriptName));
            S = [firstLine char(10) rest];
        end
