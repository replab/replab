classdef FiniteGroupType < replab.FiniteGroupType

    methods

        function iso = constructIsomorphism(self, varargin)
        % Constructs a generic isomorphism for the given elements
        %
        % Each of the argument in the variable argument list must be either:
        %
        % * a finite group whose type is equal to this finite group type,
        % * an element of this finite group type.
        %
        % Returns:
        %   `.FiniteIsomorphism`: An order preserving isomorphism whose source contains all given elements
            error('Abstract');
        end

        function G = groupWithGenerators(self, generators, varargin)
            iso = self.constructIsomorphism(generators{:});
            assert(isempty(varargin)); % TODO: support additional arguments
            G = replab.gen.FiniteGroup(generators, self, iso);
        end

    end

end
