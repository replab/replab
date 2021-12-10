classdef NiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        nElements % (integer): Number of elements in the source group
        permSet % (`+replab.+perm.PermSet`): Corresponding permutations
    end

    properties (Access = protected)
        cycloSet % (`+replab.+numerical.CycloSet`): List of matrices in the order of enumeration
        sortedToEnum % (integer(1,\*)): Permutation from sorted to enumeration convention
        enumToSorted % (integer(1,\*)): Permutation from enumeration to sorted convention
    end

    methods

        function self = NiceIsomorphism(d, sourceGenerators)
            sourceType = replab.matrix.FiniteGroupType(d);
            d = sourceType.d;
            sourceGenerators = cellfun(@(g) replab.cyclotomic(g), sourceGenerators, 'uniform', 0);
            id = replab.cyclotomic.eye(d);
            cycloSet = replab.matrix.dimino(sourceGenerators, id);
            n = cycloSet.nElements;
            targetType = replab.perm.PermutationGroupType.make(n);
            sortedToEnum = sourceType.sort(cycloSet.elements);
            enumToSorted = targetType.inverse(sortedToEnum);
            permSet = replab.perm.Set(n, n);
            perms = zeros(n, n);
            for i = 1:n
                for j = 1:n
                    perms(i,j) = enumToSorted(cycloSet.find(cycloSet.at(sortedToEnum(i)) * cycloSet.at(sortedToEnum(j))));
                end
            end
            permSet.insert(perms');
            self.permSet = permSet;
            self.cycloSet = cycloSet;
            self.sortedToEnum = sortedToEnum;
            self.enumToSorted = enumToSorted;
            self.torusMap = [];
            self.nElements = n;
            targetGenerators = cellfun(@(g) self.imageElement(g), sourceGenerators, 'uniform', 0);
            nice = replab.PermutationGroup(n, targetGenerators, 'order', n);
            self.target = nice;
            self.source = replab.MatrixGroup(d, sourceGenerators, 'nice', nice, 'niceIsomorphism', self);
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            s = replab.cyclotomic(s);
            t = self.permSet.at(self.enumToSorted(self.cycloSet.find(s)))';
        end

        % Isomorphism

        function s = preimageElement(self, t)
            s = self.cycloSet.at(self.sortedToEnum(self.permSet.find(t')));
        end

        % NiceIsomorphism

        function l = sourceContains(self, s)
            s = replab.cyclotomic(s);
            l = self.cycloSet.find(s) > 0;
        end

    end

end
