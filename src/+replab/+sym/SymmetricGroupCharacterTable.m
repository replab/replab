classdef SymmetricGroupCharacterTable < replab.CharacterTable

    properties (SetAccess = protected)
        partitions % (`.IntegerPartitions`): Integer partitions
    end

    methods

        function self = SymmetricGroupCharacterTable(n)
            partitions = replab.sym.IntegerPartitions.make(n);
            permGroup = replab.S(n);
            irrepNames = cellfun(@(part) strrep(replab.shortStr(part.conjugate), ' ', ''), partitions.list, 'uniform', 0);
            classes = partitions.conjugacyClasses;
            classNames = partitions.conjugacyClassNames;
            irreps = cellfun(@(part) permGroup.irrep(part.partition, 'specht'), partitions.list, 'uniform', 0);
            chars = zeros(length(irreps), length(classes));
            for r = 1:length(irreps)
                for c = 1:length(classes)
                    chars(r, c) = trace(irreps{r}.image(classes{c}.representative));
                end
            end
            chars = replab.cyclotomic.fromDoubles(round(chars));
            self@replab.CharacterTable(permGroup, replab.ConjugacyClasses(permGroup, classes), chars, 'irrepNames', irrepNames, 'classNames', classNames, 'irreps', irreps);

        end

    end

end
