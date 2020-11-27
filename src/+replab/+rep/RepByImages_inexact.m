classdef RepByImages_inexact < replab.RepByImages
% A finite dimensional representation of a finite group

    properties (Access = protected)
        preimages_internal
        images_internal
        imagesErrorBound_internal
    end

    methods

        function self = RepByImages_inexact(group, field, dimension, preimages, images, imagesErrorBound, varargin)
        % Constructs a representation from images of group generators and their inverses
        %
        % Keywords arguments are passed to the `+replab.Rep` constructor.
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,n) of ``group`` elements): Preimages
        %   images (cell(1,n) of double/sparse double/cyclotomic(\*,\*)): Images of the preimages
        %   imagesErrorBound (double or double(1,n) or ``[]`): Error bound on the given images
            args = struct('isUnitary', false);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            n = length(preimages);
            % if no error bound provided
            if isempty(imagesErrorBound)
                imagesErrorBound = NaN(1, n);
            end
            % if given a scalar, extend
            if isscalar(imagesErrorBound)
                imagesErrorBound = ones(1, n) * imagesErrorBound;
            end
            if any(isnan(imagesErrorBound))
                warning('Error bound missing for some/all images, computing an estimate');
            end
            % now check all images
            for i = 1:n
                img = images{i};
                if isa(img, 'replab.cyclotomic')
                    img = double(img);
                    imagesErrorBound(i) = norm(eps(img), 'fro');
                    images{i} = img;
                elseif isnan(imagesErrorBound(i))
                    eo = group.elementOrder(preimages{i});
                    imagesErrorBound(i) = norm(img^eo - eye(dimension), 'fro')/eo;
                end
            end
            self@replab.RepByImages(group, field, dimension, preimages, images, imagesErrorBound, 'isUnitary', logical(args.isUnitary), restArgs{:});
            mask = cellfun(@(g) ~group.isIdentity(g), preimages);
            self.preimages_internal = preimages(mask);
            self.images_internal = images(mask);
            self.imagesErrorBound_internal = imagesErrorBound(mask);
        end

        function c = factorization(self)
        % Returns a factorization object
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: The factorization object
            c = self.cached('factorization', @() replab.mrp.Factorization.make(self.group, self.preimages_internal, false));
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false;
        end

    end

    methods (Access = protected)

        function rho = image_double_sparse(self, g)
            rho = speye(self.dimension);
            word = self.factorization.factorize(g);
            for i = 1:length(word)
                rho = rho * self.images_internal{word(i)};
            end
        end

        function e = computeErrorBound(self)
            l = self.factorization.maximumWordLength;
            ieb = max(self.imagesErrorBound); % TODO: do better
                                              % estimate an upper bound on the operator norm (norm(X, 2))
            if self.isUnitary
                s = 1;
            else
                s = max(cellfun(@replab.numerical.norm2UpperBound, self.images_internal));
            end
            % we have (R1 + E1)*(R2 + E2)* ...*(Rl + El)
            % where Ri is the exact generator image and Ei is the error
            % now, the dominating error terms are of the form R1 R2 ... Rl-1 El
            % for some permutation of indices. We bound the error of the result using
            % norm(R1, 2) * ... * norm(Rl-1, 2) * norm(En, 'fro') up to permutation
            % which is s * ... * s * ieb, where norm(Ri, 2) <= s. There are "l" such
            % terms.
            e = s^(l-1)*ieb*l;
        end

    end

end
