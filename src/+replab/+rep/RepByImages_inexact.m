classdef RepByImages_inexact < replab.RepByImages
% A finite dimensional representation of a finite group

    properties (SetAccess = protected)
        imagesErrorBound % (double): Error bound on the given images
    end

    methods

        function self = RepByImages_inexact(group, field, dimension, preimages, images, errorBound)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/sparse double/intval/cyclotomic(\*,\*)): Images of the preimages
        %   errorBound (double): Error bound on the given images, computed as `.errorBound` on the set of images
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.isExact = false;
            % replab.Rep mutable
            n = length(preimages);
            ind = replab.mrp.inverseIndices(group, preimages);
            knownUnitary = arrayfun(@(i) ind(i) > 0 && full(all(all(images{i} == images{ind(i)}'))), 1:n);
            if knownUnitary
                self.isUnitary = true;
            end
            self.preimages = preimages;
            self.images_internal = images;
            self.imagesErrorBound = errorBound;
        end

        function c = factorization(self)
        % Returns a factorization object
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: The factorization object
            c = self.cached('factorization', @() replab.mrp.Factorization.make(self.group, self.preimages, false));
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = image_internal(self, g)
            word = self.factorization.factorize(g);
            rho = speye(self.dimension);
            for i = 1:length(word)
                rho = rho * self.images_internal{word(i)};
            end
        end

        function e = computeErrorBound(self)
            l = self.factorization.maximumWordLength;
            ieb = self.imagesErrorBound;
            % estimate an upper bound on the operator norm (norm(X, 2))
            if self.isUnitary
                s = 1;
            else
                s = max(cellfun(@replab.numerical.norm2UpperBound, self.images_internal));
            end
            % we have (R1 + E1)*(R2 + E2)* ...*(Rn + En)
            % where Ri is the exact generator image and Ei is the error
            % now, the dominating error term is of the form R1 R2 ... Rn-1 En
            % for some permutation of indices. We bound the error of the result using
            % norm(R1, 2) * ... * norm(Rn-1, 2) * norm(En, 'fro') up to permutation
            % which is s * ... * s * ieb, where norm(Ri, 2) <= s.
            e = s^(l-1)*ieb;
        end

    end

end
