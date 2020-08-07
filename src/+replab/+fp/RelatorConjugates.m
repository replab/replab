classdef RelatorConjugates < replab.Str

    properties (SetAccess = protected)
        nGenerators % (integer): Number of generators
        groupedRelators % (cell(1,\*) of cell(1,\*) of integer(1,\*)): Relators, their inverses and conjugates, grouped by first index
    end

    methods

        function self = RelatorConjugates(nGenerators, groupedRelators)
            self.nGenerators = nGenerators;
            self.groupedRelators = groupedRelators;
        end

        function res = merged(self, rhs)
            nG = self.nGenerators;
            rc = cell(1, nG*2);
            for i = 1:nG*2
                rc{i} = replab.fp.Letters.unique(horzcat(self.groupedRelators{i}, rhs.groupedRelators{i}));
            end
            res = replab.fp.RelatorConjugates(nG, rc);
        end

        function res = updated(self, newRelators)
            res = self.merged(replab.fp.RelatorConjugates.fromRelators(self.nGenerators, newRelators));
        end

    end

    methods (Static)

        function rc = fromRelators(nGenerators, relators)
            relators = cellfun(@(r) replab.fp.Letters.cyclicallyReduce(r), relators, 'uniform', 0); % line 1
            relators = relators(~cellfun(@isempty, relators));
            relators = replab.fp.Letters.unique(relators);
            relatorsC = cell(1, 0); % line 2
            for i = 1:length(relators)
                relatorsC = horzcat(relatorsC, replab.fp.Letters.cyclicConjugates(relators{i}));
                relatorsC = horzcat(relatorsC, replab.fp.Letters.cyclicConjugates(-fliplr(relators{i})));
            end
            relatorsC1 = replab.fp.Letters.unique(relatorsC);
            groupedRelators = cell(1, nGenerators*2);
            for i = 1:nGenerators*2
                groupedRelators{i} = cell(1, 0);
            end
            for i = 1:length(relatorsC1)
                r = relatorsC1{i};
                f = r(1);
                if f < 0
                    f = nGenerators + (-f);
                end
                groupedRelators{f} = horzcat(groupedRelators{f}, r);
            end
            rc = replab.fp.RelatorConjugates(nGenerators, groupedRelators);
        end

    end

end
