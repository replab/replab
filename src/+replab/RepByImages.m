classdef RepByImages < replab.Rep
% A finite dimensional representation of a finite group

    properties (SetAccess = protected)
        preimages % (cell(1,\*) of `.group` elements): Preimages
        images % (cell(1,\*) of double(\*,\*) or `.cyclotomic`, may be sparse): Images
    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'preimages';
            names{1, end+1} = 'images';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:length(self.images)
                names{1, end+1} = sprintf('preimages{%d}', i);
                values{1, end+1} = self.preimages{i};
                names{1, end+1} = sprintf('images{%d}', i);
                values{1, end+1} = self.images{i};
            end
        end

        % Rep


        function res = imap(self, f)
            preimages1 = cellfun(@(p) f.imageElement(p), self.preimages, 'uniform', 0);
            res = replab.RepByImages.make(f.image, self.field, self.dimension, preimages1, self.images);
            res.isExact = self.isExact;
            res.isUnitary = self.isUnitary;
            res.trivialDimension = self.trivialDimension;
            res.isIrreducible = self.isIrreducible;
            res.frobeniusSchurIndicator = self.frobeniusSchurIndicator;
            res.isDivisionAlgebraCanonical = self.isDivisionAlgebraCanonical;
        end

    end

    methods (Static)

        function rep = make(group, field, dimension, preimages, images)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double(\*,\*), may be sparse or symbolic): Images of the preimages
            rep = replab.rep.repByImages(group, field, dimension, preimages, images);
        end

        function rep1 = fromExactRep(rep)
        % Constructs a `.RepByImages` from an existing representation with exact images
            assert(isa(rep.group, 'replab.NiceFiniteGroup'));
            nG = rep.group.nGenerators;
            images = arrayfun(@(i) rep.image_internal(rep.group.generator(i)), 1:nG, 'uniform', 0);
            rep1 = replab.RepByImages.make(rep.group, rep.field, rep.dimension, rep.group.generators, images);
        end

        function rep = fromImageFunction(group, field, dimension, imageFun)
        % Constructs a RepByImages representation using a given morphism
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group of which to construct a representation
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Dimension of the representation
        %   imageFun (function_handle): Function that returns a matrix for any element of ``G``
        %
        % Returns:
        %   `+replab.RepByImages`: The constructed representation
            assert(isa(group, 'replab.FiniteGroup'), 'The given group must be a FiniteGroup');
            nG = group.nGenerators;
            images = arrayfun(@(i) imageFun(group.generator(i)), 1:nG, 'uniform', 0);
            rep = replab.RepByImages.make(group, field, dimension, group.generators, images);
        end

    end

end
