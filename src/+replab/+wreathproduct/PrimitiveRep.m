classdef PrimitiveRep < replab.Rep
% Primitive representation of a wreath product group
%
% See `replab.wreathproduct.Common.primitiveRep`
%
% See `replab.wreathproduct.ImprimitiveRep`
    
    properties (SetAccess = protected)
        Arep % replab.Rep: Representation of the group whose copies are acted upon
    end
    
    methods
        
        function self = PrimitiveRep(group, Arep)
            assert(isa(group, 'replab.wreathproduct.Common'));
            assert(group.A == Arep.group);
            dA = Arep.dimension;
            n = group.H.domainSize;
            self.Arep = Arep;
            self.dimension = dA^n;
            self.isUnitary = Arep.isUnitary;
            self.field = Arep.field;
            self.group = group;
        end
        
        function rho = image(self, g)
            h = g{1};
            base = g{2};
            rho = self.Arep.image(base{1});
            n = self.group.H.domainSize;
            dA = self.Arep.dimension;
            for i = 2:n
                rho = kron(rho, self.Arep.image(base{i}));
            end
            dims = dA * ones(1, n);
            d = self.dimension;
            I = permute(reshape(1:d, dims), fliplr(n + 1 - h));
            I = I(:)';
            rho = sparse(I, 1:d, ones(1, d), d, d) * rho;
            if ~replab.Settings.useSparse
                rho = full(rho);
            end
        end
        
    end
    
end
