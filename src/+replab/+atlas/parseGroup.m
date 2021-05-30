function group = parseGroup(data)
% Parses group information from the data from a JSON file
%
% Valid JSON looks like::
%
%  {
%    relators: [
%      "x1^3",
%      "x2^3",
%      "x1 x2 x1 x2"
%    ],
%    order: "12",
%    generatorNames: ["x1", "x2"],
%    permutationGenerators: [
%      [4, 1, 3, 2],
%      [4, 2, 1, 3]
%    ],
%    classes: [
%      [1, 2, 3, 4],
%      [2, 1, 4, 3],
%      [3, 1, 2, 4],
%      [4, 1, 3, 2]
%    ],
%    derivedSeriesAbelianInvariants: [[3], [2, 2], []]
%  }
%
% In the parsed data, arrays are always cell arrays.
%
% Note that relators can be given as integer arrays as well (thus ``[[1 1 1], [2 2 2], [1 2 1 2]]`` for the example above).
%
%
% Example:
%   >>> json = '{"relators" : [[1,1,1],[2,2,2],[1,2,1,2]],"order" : "12","generatorNames" : ["x1","x2"],"permutationGenerators" : [[4,1,3,2],[4,2,1,3]],"classes" : [[1,2,3,4],[2,1,4,3],[3,1,2,4],[4,1,3,2]],"derivedSeriesAbelianInvariants" : [[3], [2,2], []]}';
%   >>> data = replab.util.parseJSON(json);
%   >>> G = replab.atlas.parseGroup(data);
%
% Args:
%   data (struct): Struct read from a JSON file
%
% Returns:
%   `.AbstractGroup`: Abstract group data
    generatorNames = data.generatorNames;
    relators = data.relators;
    for i = 1:length(relators)
        if isa(relators{i}, 'cell')
            relators{i} = replab.fp.Letters.print(cell2mat(relators{i}), generatorNames);
        end
    end
    order = data.order;
    switch class(order)
      case 'double'
        assert(order <= 2^53, 'Write orders > 2^53 as strings');
        order = vpi(order);
      case 'char'
        order = vpi(order);
      otherwise
        error('Incorrect order type');
    end
    permutationGenerators = cellfun(@cell2mat, data.permutationGenerators, 'uniform', 0);
    group = replab.AbstractGroup(generatorNames, relators, 'permutationGenerators', permutationGenerators);
    nClasses = length(data.classes);
    classes = cell(1, nClasses);
    for i = 1:nClasses
        cd = data.classes{i};
        if isa(cd, 'char')
            classes{i} = replab.ConjugacyClass.make(group, cd);
        else
            classes{i} = replab.ConjugacyClass.make(group, group.niceMorphism.preimageElement(cell2mat(cd)));
        end
    end
    group.cache('conjugacyClasses', replab.ConjugacyClasses(group, classes), 'error');
end
