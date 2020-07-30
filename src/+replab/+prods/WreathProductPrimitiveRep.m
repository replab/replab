classdef WreathProductPrimitiveRep < replab.Rep
% Primitive representation of a wreath product group

    properties (SetAccess = protected)
        Arep % (`+replab.Rep`): Representation of the group whose copies are acted upon
    end

    methods

        function self = WreathProductPrimitiveRep(group, Arep)
        % Constructs a primitive representation of a wreath product group
        %
        % Args:
        %   group (`+replab.WreathProductGroup`): Wreath product group
        %   Arep (`+replab.Rep`): Representation of the wreath product base factor
        %
        % Returns:
        %   `+replab.Rep`: A wreath product group representation
            assert(isa(group, 'replab.WreathProductGroup'));
            assert(group.A == Arep.group);
            dA = Arep.dimension;
            n = group.H.domainSize;
            self.Arep = Arep;
            self.dimension = dA^n;
            self.isUnitary = Arep.isUnitary;
            self.field = Arep.field;
            self.group = group;
        end

        function rho = image_internal(self, g)
            h = g{1};
            base = g{2};
            rho = self.Arep.image_internal(base{1});
            n = self.group.H.domainSize;
            dA = self.Arep.dimension;
            for i = 2:n
                rho = kron(rho, self.Arep.image_internal(base{i}));
            end
            dims = dA * ones(1, n);
            d = self.dimension;
            I = permute(reshape(1:d, dims), fliplr(n + 1 - h));
            I = I(:)';
            rho = sparse(I, 1:d, ones(1, d), d, d) * rho;
            rho = full(rho);
        end

    end

end
