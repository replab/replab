classdef FiniteGroupType < replab.FiniteGroupType

    methods (Access = protected)

        function G = makeGenericGroup(self, generators, nice, niceIsomorphism)
        % Creates a generic group of this type
        %
        % Args:
        %   generators (cell(1,\*) of elements of this type): Group generators
        %   nice (`+replab.FiniteGroup`): Subgroup of the isomorphism target whose generators are in 1-to-1 correspondance with ``generators``
        %   niceIsomorphism (`+replab.+gen.NiceIsomorphism`): Isomorphism whose source contains all the generators
            G = replab.gen.FiniteGroup(self, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
        end
    end

    methods

        function iso = constructIsomorphism(self, elements)
        % Constructs a generic isomorphism for the given elements
        %
        % Args:
        %   elements (cell(1,\*) of group type elements): Elements
        %
        % Returns:
        %   `.FiniteIsomorphism`: An order preserving isomorphism whose source contains all given elements
            error('Abstract');
        end

    end

    methods % Implementations

        function G = groupWithGenerators(self, generators, varargin)
            assert(all(ismember(varargin(1:2:end), {'generatorNames', 'order', 'relators'})));
            niceIso = self.constructIsomorphism(generators);
            target = niceIso.target;
            targetGenerators = cellfun(@(g) niceIso.imageElement(g), generators, 'uniform', 0);
            nice = target.subgroup(targetGenerators, varargin{:});
            G = self.makeGenericGroup(generators, nice, niceIso);
        end

    end

end
