classdef RepByImagesFactorization < replab.Rep
% Constructs a representation from images

    properties
        preimages % (cell(1,n) of `.group` elements): Preimages, the set is closed under the group inverse operation
        factorization % (`+replab.+mrp.Factorization`): Factorization of group elements
        approxImages % (cell(1,n) of double(d,d)): Approximate images
        errorImages2 % (double(1,n)): Upper bound on the image error in 2-norm
        errorImagesF % (double(1,n)): Upper bound on the image error in Frobenius norm
        normImages2 % (double(1,n)): Upper bound on the image 2-norm
    end

    methods (Static)

        function ok = areImagesUnitary(group, preimages, images)
        % Returns whether the approximative images define a unitary representation
        %
        % By convention, this is the case when the set ``preimages`` is closed under the group
        % inverse operation, and when ``images{i}`` is the conjugate transpose of ``images{inv(i)}``, defining
        % ``inv(i)`` such that ``preimages{i}`` is the inverse of ``preimages{inv(i)}``.
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Group to check the images of
        %   preimages (cell(1,n) of elements of ``group``): Preimages
        %   images (cell(1,n) of double(d,d)): Approximate images
            inv = replab.mrp.inverseIndices(group, preimages);
            for i = 1:length(preimages)
                if i <= inv(i)
                    img_i = images{i};
                    img_inv = images{inv(i)};
                    if any(any(img_i ~= img_inv'))
                        ok = false;
                        return
                    end
                end
            end
            ok = true;
        end

    end

    methods

        function self = RepByImagesFactorization(group, field, dimension, preimages, approxImages, errorImages2, errorImagesF)
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            n = length(self.preimages);
            self.preimages = preimages;
            self.approxImages = approxImages;
            self.errorImages2 = errorImages2;
            self.errorImagesF = errorImagesF;
            self.factorization = replab.mrp.Factorization.make(group, preimages, false);
            if replab.rep.RepByImagesFactorization.areImagesUnitary(group, preimages, approxImages)
                self.isUnitary = true;
                self.normImages2 = errorImages2 + 1; % 2-norm of a unitary matrix is 1
            else
                self.normImages2 = cellfun(@(X) norm(X, 2), approxImages);
            end
        end

        function rho = image_internal(self, g)
            rho = self.image_error(g);
        end

        function [rho error2 errorF] = image_error(self, g)
            w = self.factorization.preimageElement(g);
            l = w(1);
            error2 = self.errorImages2(l);
            errorF = self.errorImagesF(l);
            norm2 = self.normImages2(l);
            rho = self.approxImages{l};
            for i = 2:length(w)
                l = w(i);
                rho = rho * self.approxImages{l};
                error2 = error2 * self.normImages2(l) + self.errorImages2(l) * norm2;
                normF = norm2 * sqrt(self.dimension); % upper bound
                errorF = errorF * self.normImagesF(l) + self.errorImagesF(l) * norm2 * sqrt(self.dimension);
            end
        end

    end

end
