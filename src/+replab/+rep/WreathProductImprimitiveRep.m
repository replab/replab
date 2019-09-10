classdef WreathProductImprimitiveRep < replab.Rep
    
    properties (SetAccess = protected)
        Arep;
    end
    
    methods
        
        function self = WreathProductImprimitiveRep(group, Arep)
            assert(isa(group, 'replab.WreathProductGroup'));
            assert(isequal(group.A, Arep.group));
            dA = Arep.dimension;
            n = group.H.domainSize;
            self.Arep = Arep;
            self.dimension = n*dA;
            self.field = Arep.field;
            self.group = group;
        end
        
        function rho = image(self, g)
            h = g{1};
            base = g{2};
            n = length(h);
            rhos = arrayfun(@(i) self.Arep.image(base{i}), 1:n, 'uniform', 0);
            rho = blkdiag(rhos{:});
            dA = self.Arep.dimension;
            rho = kron(sparse(h, 1:n, ones(1,n), n, n), speye(dA)) * rho;
            if ~replab.Settings.useSparse
                rho = full(rho);
            end
        end
        
    end
    
end
