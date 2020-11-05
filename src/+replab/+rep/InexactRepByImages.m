classdef InexactRepByImages < replab.RepByImages
% A finite dimensional representation of a finite group

    methods

        function self = InexactRepByImages(group, field, dimension, preimages, images, knownUnitary)
        % Constructs a representation from images of group generators and their inverses
        %
        % Args:
        %   group (`+replab.FiniteGroup`): Finite group represented
        %   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   preimages (cell(1,\*) of ``group`` elements): Preimages
        %   images (cell(1,\*) of double/sparse double/intval/cyclotomic(\*,\*)): Images of the preimages
        %   knownUnitary (logical, optional): Whether the representation is known, defaults to false
            if nargin < 6 || isempty(knownUnitary)
                knownUnitary = false;
            end
            % replab.Rep immutable
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            % replab.Rep mutable
            if knownUnitary
                self.isUnitary = true;
            end
            self.preimages = preimages;
            self.images = images;
        end

        function c = factorization(self)
        % Returns a factorization object
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: The factorization object
            c = self.cached('factorization', @() replab.mrp.Factorization.make(self.group, self.preimages, false))
        end

    end

    methods % Implementations

        % Rep

        function rho = image_internal(self, g)
            word = self.factorization.factorize(g);
            rho = speye(self.dimension);
            for i = 1:length(word)
                rho = rho * self.images{word(i)};
            end
        end

    end

end
