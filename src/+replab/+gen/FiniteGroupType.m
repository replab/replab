classdef FiniteGroupType < replab.FiniteGroupType

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

        function G = groupWithGenerators(self, generators, varargin)
            assert(all(ismember(varargin(1:2:end), {'generatorNames', 'order', 'relators'})));
            iso = self.constructIsomorphism(generators);
            targetGenerators = cellfun(@(g) iso.imageElement(g), generators, 'uniform', 0);
            nice = iso.target.subgroup(targetGenerators, varargin{:});
            G = replab.gen.FiniteGroup(self, nice, iso, generators);
        end

    end

end
