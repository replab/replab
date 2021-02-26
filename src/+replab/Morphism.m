classdef Morphism < replab.Obj
% Describes a morphism between groups
%
% If both groups are compact groups with a `+replab.CompactGroup.reconstruction` available, the `.torusMap` linear map
% expresses the relationship between elements of the torus of the source and the target. More specifically, let ``s``
% be an element of the torus ``.source.reconstruction.source``, there is a corresponding element ``t`` of the torus
% ``.target.reconstruction.source`` such that ``t = s * torusMap``.
    properties (SetAccess = protected)
        source % (`.Group`): Source group
        target % (`.Group`): Target group
        torusMap % (integer(\*,\*) or ``[]``): Relation between the tori of source and target if known
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

        function t = imageTorusElement(self, s)
        % Returns the image of an element of the source torus (used in optimizations)
        %
        % Args:
        %   s (element of ``source.reconstruction.source``): Torus element of source
        %
        % Returns:
        %   element of ``target.reconstruction.source``): Torus element of target
            assert(self.source.hasReconstruction && self.target.hasReconstruction);
            t = mod(s * self.torusMap, 1);
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

        function m = lambda(source, target, imageElementFun, torusMap)
        % Creates a morphism from an image function
        %
        % Args:
        %   source (`.Group`): Source group
        %   target (`.Group`): Target group
        %   imageElementFun (function_handle): Function computing images of elements
        %   torusMap (integer(\*,\*)): Torus map to use in the morphism construction
        %
        % Returns:
        %   `+replab.Morphism`: Constructed morphism
            if nargin < 4
                torusMap = [];
            end
            m = replab.mrp.Lambda(source, target, imageElementFun, torusMap);
        end

    end

end
