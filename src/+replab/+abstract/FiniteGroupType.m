classdef FiniteGroupType < replab.gen.StaticFiniteGroupType

    properties (SetAccess = protected)
        generatorNames % (cell(1,\*) of charstring): Generator names
        typeId % (integer): Unique abstract group type ID
        permutationGroup % (`+replab.PermutationGroup`): Concrete instance of the abstract group
    end

    methods

        function self = FiniteGroupType(generatorNames, permutationGroup)
            self.generatorNames = generatorNames;
            self.typeId = replab.globals.nextUniqueId;
            self.permutationGroup = permutationGroup;
            self.niceType = permutationGroup.type;
            self.isomorphism = replab.gen.StaticNiceIsomorphism(self, generatorNames, permutationGroup);
            self.identity = '1';
            % we skip finishConstruction because the target group is already known
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(self.imageElement(x) == self.imageElement(y));
        end

        function s = sample(self)
            t = self.isomorphism.target.sample;
            s = self.preimageElement(t);
        end

        % Monoid

        function z = compose(self, x, y)
            xl = replab.fp.Letters.parse(x, self.generatorNames);
            yl = replab.fp.Letters.parse(y, self.generatorNames);
            zl = replab.fp.Letters.compose(xl, yl);
            z = replab.fp.Letters.print(zl, self.generatorNames);
        end

        % Group

        function z = inverse(self, x)
            xl = replab.fp.Letters.parse(x, self.generatorNames);
            zl = replab.fp.Letters.inverse(xl);
            z = replab.fp.Letters.print(zl, self.generatorNames);
        end

        % FiniteGroupType

        function l = isSameTypeAs(self, otherType)
            l = isa(otherType, 'replab.abstract.FiniteGroupType') && self.typeId == otherType.typeId;
        end

        % StaticFiniteGroupType

        function t = imageElement(self, s)
            t = self.permutationGroup.imageWord(s);
        end

        function S = makeParentGroup(self, generators, nice, niceIsomorphism)
            S = []; % will be filled up afterwards
        end

        function s = preimageElement(self, t)
            s = self.permutationGroup.factorizeWord(t);
        end

    end

end
