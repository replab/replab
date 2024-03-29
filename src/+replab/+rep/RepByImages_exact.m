classdef RepByImages_exact < replab.RepByImages
% A finite dimensional representation of a finite group
%
% The finite group must have an isomorphism to a permutation group, so that the images of that representation
% can be reprseented using a BSGS construction that storing the stabilizer chain with transversal elements both encoding the
% group transversals and their images (see `+replab.+bsgs.ChainWithImages`).

    properties (SetAccess = protected)
        inverseImages % (cell(1,n) of double(d,d) or cyclotomic(d,d), may be sparse): Images
    end

    methods

        function self = RepByImages_exact(group, field, dimension, preimages, images, varargin)
        % Constructs a representation from images of group generators and their inverses
        %
        % Keywords arguments are passed to the `+replab.Rep` constructor.
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/sparse double/cyclotomic(\*,\*)): Images of the preimages
            inverseImages = replab.rep.computeInverses(group, @(x,y) x*y, preimages, images);
            imagesErrorBound = zeros(1, length(preimages));
            isUnitary = all(arrayfun(@(i) full(all(all(images{i} == inverseImages{i}'))), 1:length(preimages)));
            [args, exists, oldValue] = replab.util.keyValuePairsUpdate(varargin, 'isUnitary', isUnitary);
            assert(~exists || isempty(oldValue) || isequal(oldValue, isUnitary), 'Wrong isUnitary keyword argument value');
            self@replab.RepByImages(group, field, dimension, preimages, images, imagesErrorBound, args{:});
            self.inverseImages = inverseImages;
        end

        function c = chain(self)
        % Returns a BSGS chain that computes images of this representation
        %
        % The preimages are permutations.
        %
        % Returns:
        %   `+replab.+bsgs.ChainWithImages`: The BSGS chain with cyclotomic matrix images
            c = self.cached('chain', @() self.computeChain);
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = true;
        end

    end

    methods (Access = protected)

        function res = computeChain(self)
            preimages = self.preimages;
            images = self.images;
            inverseImages = self.inverseImages;
            isId = cellfun(@(g) self.group.isIdentity(g), preimages);
            preimages = preimages(~isId);
            images = images(~isId);
            inverseImages = inverseImages(~isId);
            m = length(preimages);
            d = self.dimension;
            iso = self.group.permutationIsomorphism;
            n = self.group.permutationIsomorphism.target.domainSize;
            nicePreimages = cellfun(@(g) iso.imageElement(g), self.preimages, 'uniform', 0);
            order = self.group.order;
            base = [];
            if self.isUnitary
                % for unitary/orthogonal representations, we don't need to compute inverses explicitly
                if self.overR
                    target = replab.O(d);
                else
                    target = replab.U(d);
                end
                res = replab.bsgs.ChainWithImages.make(n, target, nicePreimages, images, base, order);
            else
                % for nonunitary/nonorthogonal representations, we store the inverses alongside the matrix, and use a special
                % group structure to perform computations
                target1 = replab.domain.GeneralLinearGroupWithInverses(self.field, self.dimension, true);
                target2 = replab.domain.GeneralLinearGroup(self.field, self.dimension, true);
                cut = replab.Morphism.lambda(target1, target2, @(X) X(:, 1:self.dimension));
                images = arrayfun(@(i) [images{i} inverseImages{i}], 1:m, 'uniform', 0);
                res = replab.bsgs.ChainWithImages.make(n, target1, nicePreimages, images, base, order);
                res = res.mapImages(cut);
            end
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = image_double_sparse(self, g)
            perm = self.group.permutationIsomorphism.imageElement(g);
            rho = self.chain.image(perm);
            if isa(rho, 'replab.cyclotomic')
                rho = double(rho);
            end
        end

        function rho = image_exact(self, g)
            perm = self.group.permutationIsomorphism.imageElement(g);
            rho = self.chain.image(perm);
            if isa(rho, 'double')
                rho = replab.cyclotomic(rho);
            end
        end

        function e = computeErrorBound(self)
            if self.isUnitary
                % coefficients are bounded in [-1,1], so their error is at most eps(1)
                % and there is dim*dim of them, under the square root
                e = eps(1)*self.dimension;
                if self.overC
                    e = e*sqrt(2); % real and imaginary part
                end
            else
                simRep = self.unitarize;
                % worst case is given by the condition number of the change of basis
                e = eps(1)*self.dimension*self.conditionNumberEstimate;
                if self.overC
                    e = e*sqrt(2); % real and imaginary part
                end
            end
        end

    end

end
