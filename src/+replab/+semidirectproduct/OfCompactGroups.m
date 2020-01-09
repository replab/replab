classdef OfCompactGroups < replab.CompactGroup
% Describes an external semidirect product of compact groups
%
% See `+replab.CompactGroup.semidirectProduct`
    
    properties (SetAccess = protected)
        H % replab.CompactGroup: Group acting
        N % replab.CompactGroup: Group acted upon
        phi % replab.Action: Action of H on N
    end

    methods
        
        function self = OfCompactGroups(phi)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            self.phi = phi;
            assert(isa(H, self.requiredType));
            assert(isa(N, self.requiredType));
            self.H = H;
            self.N = N;
            self.identity = {H.identity N.identity};
        end
        
        function t = requiredType(self)
            t = 'replab.CompactGroup';
        end
        
        %% Domain methods
        
        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end
        
        function g = sample(self)
            g = {self.H.sample self.N.sample};
        end
        
        %% Monoid methods

        function z = compose(self, x, y)
        % Composition
        %
        % Relation to phi is the conjugation
        % phi_h(n) = h n h^-1
        % we have z = xh xn yh yn = xh yh yh^-1 xn yh yn =
        % = xh yh phi_(yh^-1)(xn) yn
        % and thus
        % zh = xh yh
        % zn = phi_(yh^-1)(xn) yn
            xh = x{1};
            xn = x{2};
            yh = y{1};
            yn = y{2};
            yhinv = self.H.inverse(yh);
            zh = self.H.compose(xh, yh);
            zn = self.N.compose(self.phi.leftAction(yhinv, xn), yn);
            z = {zh zn};
        end
        
        %% Group methods
        
        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = self.N.inverse(self.phi.leftAction(xh, xn));
            z = {zh zn};
        end
        
        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = {self.H.sampleUniformly self.N.sampleUniformly};
        end
        
    end

end
