classdef NiceIsomorphism < replab.FiniteIsomorphism
% An isomorphism from a subgroup of a finite group type to a finite group (Abstract)
%
% In addition, the isomorphism should preserve the type order.
%
% Its source must be of a generic type `+replab.+gen.FiniteGroupType`.

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.NiceIsomorphismLaws(self);
        end

        % FiniteIsomorphism

        function S = preimageGroup(self, T)
            S = self.source.type.groupFromNiceImage(T, self);
        end

        function l = preservesTypeOrder(self)
            l = true;
        end

        function T = imageGroup(self, S)
            if S.compatibleWithNiceIsomorphism(self)
                T = S.nice;
            else
                T = imageGroup@replab.FiniteIsomorphism(self, S);
            end
        end

    end

    methods

        function l = sourceContains(self, s)
        % Returns whether the source of this morphism contains the given type element
        %
        % Args:
        %   s (element of ``.source.type``): Element to check
        %
        % Returns:
        %   logical: True if source contains the given element
            error('Abstract');
        end

    end

end
