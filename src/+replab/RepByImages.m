classdef RepByImages < replab.Rep
% A finite dimensional representation of a finite group

    properties (SetAccess = protected)
        preimages % (cell(1,n) of `.group` elements): Preimages
        images % (cell(1,n) of double(d,d) or cyclotomic(d,d) or intval(d,d), may be sparse): Images
        imagesErrorBound % (double(1,n)): Error bound on the given images
    end

    methods

        function self = RepByImages(group, field, dimension, preimages, images, imagesErrorBound, varargin)
            self@replab.Rep(group, field, dimension, varargin{:});
            self.preimages = preimages;
            self.images = images;
            self.imagesErrorBound = imagesErrorBound;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_isTrivial(self, options)
        % Rewrite rule: replace this representation by a trivial representation if it is a representation of the trivial group
            if length(self.images) == 0
                res = self.group.trivialRep(self.field, self.dimension);
            else
                res = [];
            end
        end

        function res = rewriteTerm_hasBlockStructure(self, options)
        % Rewrite rule: rewrite this representation as a direct sum if it has a block structure
            if self.isExact && length(self.images) > 0
                n = length(self.images);
                mask = self.images{1} ~= 0;
                for i = 2:length(self.images)
                    mask = mask | (self.images{i} ~= 0);
                end
                G = replab.UndirectedGraph.fromAdjacencyMatrix(double(mask));
                P = G.connectedComponents;
                blocks = P.blocks;
                if length(blocks) == 1
                    res = [];
                    return
                end
                m = length(blocks);
                factors = cell(1, m);
                for i = 1:m
                    block = blocks{i};
                    blockImages = cell(1, n);
                    for j = 1:n
                        img = self.images{j};
                        blockImages{j} = img(block, block);
                    end
                    factors{i} = self.group.repByImages(self.field, length(block), ...
                                                        'preimages', self.preimages, 'images', blockImages);
                end
                ds = self.group.directSumRep(self.field, factors);
                perm = [blocks{:}];
                D = self.dimension;
                A = sparse(perm, 1:D, ones(1, D), D, D);
                Ainv = sparse(1:D, perm, ones(1, D), D, D);
                res = ds.similarRep(A, 'inverse', Ainv);
            else
                res = [];
            end
        end

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

% $$$         function res = imap(self, f)
% $$$             preimages1 = cellfun(@(p) f.imageElement(p), self.preimages, 'uniform', 0);
% $$$             res = replab.RepByImages.make(f.image, self.field, self.dimension, preimages1, self.images_internal);
% $$$             res.isExact = self.isExact;
% $$$             res.isUnitary = self.isUnitary;
% $$$             res.trivialDimension = self.trivialDimension;
% $$$             res.isIrreducible = self.isIrreducible;
% $$$             res.frobeniusSchurIndicator = self.frobeniusSchurIndicator;
% $$$             res.isDivisionAlgebraCanonical = self.isDivisionAlgebraCanonical;
% $$$         end

    end

    methods (Static)

        function rep1 = fromExactRep(rep)
        % Constructs a `+replab.RepByImages` from an existing representation with exact images
            assert(isa(rep.group, 'replab.FiniteGroup'));
            assert(rep.isExact);
            gens = rep.group.generators;
            images = cellfun(@(g) rep.image_internal(g), gens, 'uniform', 0);
            rep1 = rep.group.repByImages(rep.field, rep.dimension, 'preimages', gens, 'images', images);
            rep1.copyProperties(rep);
        end
% $$$
% $$$         function rep = fromImageFunction(group, field, dimension, imageFun)
% $$$         % Constructs a RepByImages representation using a given morphism
% $$$         %
% $$$         % Args:
% $$$         %   group (`+replab.FiniteGroup`): Group of which to construct a representation
% $$$         %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
% $$$         %   dimension (integer): Dimension of the representation
% $$$         %   imageFun (function_handle): Function that returns a matrix for any element of ``G``
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.RepByImages`: The constructed representation
% $$$             assert(isa(group, 'replab.FiniteGroup'), 'The given group must be a FiniteGroup');
% $$$             nG = group.nGenerators;
% $$$             images = arrayfun(@(i) imageFun(group.generator(i)), 1:nG, 'uniform', 0);
% $$$             rep = replab.RepByImages.make(group, field, dimension, group.generators, images);
% $$$         end
    end

end
