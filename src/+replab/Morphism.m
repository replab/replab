classdef Morphism < replab.Obj
% Describes a morphism between groups

    properties (SetAccess = protected)
        source % (`.Group`): Source group
        target % (`.Group`): Target group
    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.MorphismLaws(self);
        end

    end

    methods % Images

        function t = imageElement(self, s)
        % Returns the image of the given source element
        %
        % Args:
        %   s (element of `.source`): Element to compute the image of
        %
        % Returns:
        %   element of `.target`: Image
            error('Abstract');
        end

    end

    methods % Morphism composition

        function res = compose(self, applyFirst)
        % Composition of morphisms, the right hand side applied first
        %
        % Note: if both morphisms are finite morphisms, the resulting morphism is also a finite morphism.
        % Note: if both morphisms are isomorphisms, the resulting morphism is also an isomorphism.
        %
        % Args:
        %   applyFirst (`.Morphism`): Morphism to apply first
        %
        % Returns:
        %   `.Morphism`: The composition of the given morphism applied first, followed by this morphism.
            res = replab.mrp.compose(self, applyFirst);
        end

        function res = andThen(self, applyLast)
        % Composition of morphisms, the right hand side applied last
        %
        % Args:
        %   applyLast (`.Morphism`): Morphism to apply last
        %
        % Returns:
        %   `.Morphism`: The composition of this morphism applied first, followed by the given morphism
            res = replab.mrp.compose(applyLast, self);
        end

        function res = mtimes(self, applyFirst)
        % Shorthand for `.compose`
            res = self.compose(applyFirst);
        end

    end

    methods (Static) % Morphism creation

        function m = lambda(source, target, imageElementFun)
        % Creates a morphism from an image function
        %
        % Args:
        %   source (`.Group`): Source group
        %   target (`.Group`): Target group
        %   imageElementFun (function_handle): Function computing images of elements
        %
        % Returns:
        %   `+replab.Morphism`: Constructed morphism
            m = replab.mrp.Lambda(source, target, imageElementFun);
        end

    end

end
