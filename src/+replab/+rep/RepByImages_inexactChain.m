classdef RepByImages_inexactChain < replab.RepByImages
% A finite dimensional representation of a finite group
%
% Obtained by an approximation of an exact BSGS chain

    properties (SetAccess = protected)
        chain % (`+replab.+bsgs.Chain`): BSGS chain with approximated images
    end

    methods

        function self = RepByImages_inexactChain(rep)
        % Constructs a representation from images of group generators and their inverses
        %
        % Keywords arguments are passed to the `+replab.Rep` constructor.
        %
        % Args:
        %   rep (`+replab.+rep.RepByImages_exact`): Exact representation
            preimages = rep.preimages;
            n = length(preimages);
            imagesErrorBound = zeros(1, n);
            images = cell(1, n);
            for i = 1:n
                img = rep.images{i};
                if isa(img, 'replab.cyclotomic')
                    [res, error] = img.doubleApproximation;
                    images{i} = res;
                    imagesErrorBound(i) = norm(error, 'fro');
                else
                    images{i} = img;
                    imagesErrorBound(i) = 0;
                end
            end
            args = rep.knownProperties({'isUnitary', 'trivialDimension', 'frobeniusSchurIndicator', 'isDivisionAlgebraCanonical', 'kernel'});
            self@replab.RepByImages(rep.group, rep.field, rep.dimension, preimages, images, imagesErrorBound, args{:});
            [self.chain, errorBound] = rep.chain.double(rep.knownUnitary);
            self.cache('errorBound', errorBound, 'error');
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false;
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = computeDouble(self)
            rho = self;
        end

        function rho = image_double_sparse(self, g)
            perm = self.group.niceMorphism.imageElement(g);
            rho = self.chain.image(perm);
        end

    end

end
