classdef Morphism < replab.Str
% Describes a morphism between groups

    properties (SetAccess = protected)
        source % (`.Group`): Source group
        target % (`.Group`): Target group
    end

    methods

        function t = image(s)
            error('Abstract');
        end

    end

    methods (Static)

        function m = lambda(source, target, imageFun)
            m = replab.mrp.LambdaMorphism(source, target, imageFun);
        end

        function m = compose(second, first)
            m = replab.mrp.CompositionMorphism(second, first);
        end

        function m = byImages(source, target, generatorImages)
        % Constructs a morphism of a NiceFiniteGroup from images of its generators
        %
        % Args:
        %   source (`.NiceFiniteGroup`): Source of the morphism
        %   target (`.Group`): Target of the morphism
        %   generatorImages (cell(1,\*) of target elements): Images of the generators of ``source``
        %
        % Returns:
        %   `.Morphism`: The constructed morphism
            if isa(source, 'replab.PermutationGroup')
                if isa(target, 'replab.PermutationGroup')
                    m = replab.mrp.PermPermMorphism.byImages(source, generatorImages);
                    return
                else
                    m = replab.mrp.PermMorphism.byImages(source, target, generatorImages);
                    return
                end
            elseif isa(source, 'replab.NiceFiniteGroup')
                first = source.niceMonomorphism;
                second = replab.Morphism.make(source.niceGroup, target, generatorImages);
                m = replab.mrp.CompositionMorphism(second, first);
            end
            error('Unsupported');
        end

    end

end
