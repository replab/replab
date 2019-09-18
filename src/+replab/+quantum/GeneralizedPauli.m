classdef GeneralizedPauli < replab.FiniteGroup
% Returns a definition of the generalized Pauli group
%
% See for example https://arxiv.org/pdf/quant-ph/0408190.pdf
% 
% An element of the Pauli group of dimension d is represented by
% three integers a, b, c = 0,...,d-1 in a row vector [a b c]
% such that for g = [a b c], the represented group element is
%    w^a x^b z^c
%
% w is a common phase
% x is a cyclic shift of levels
% z is a level-dependent phase 

    properties
        d; % Dimension of the qudit
    end
        
    methods
        
        function self = GeneralizedPauli(d)
            self.d = d;
            self.identity = [0 0 0];
            % group generators
            w = [1 0 0];
            x = [0 1 0];
            z = [0 0 1];
            self.generators = {w x z};
        end
        
        function o = order(self)
            o = vpi(self.d)^3;
        end
        
        function E = elements(self)
            E = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.elementAt(ind), ...
                                            @(el) self.elementFind(el));
        end
        
        function x3 = compose(self, x1, x2)
            d = self.d;
            a1 = x1(1);
            b1 = x1(2);
            c1 = x1(3);
            a2 = x2(1);
            b2 = x2(2);
            c2 = x2(3);
            a3 = mod(a1 + a2 + c1*b2, d);
            b3 = mod(b1 + b2, d);
            c3 = mod(c1 + c2, d);
            x3 = [a3 b3 c3];
        end
        
        function b = eqv(self, x1, x2)
            b = isequal(x1, x2);
        end
        
        function x = sample(self)
            x = randi([0 self.d-1], 1, 3);
        end
        
        function x2 = inverse(self, x1)
            d = self.d;
            a1 = x1(1);
            b1 = x1(2);
            c1 = x1(3);
            a2 = mod(d - a1 + b1*c1, d);
            b2 = mod(d - b1, d);
            c2 = mod(d - c1, d);
            x2 = [a2 b2 c2];
        end
        
        function D = decomposition(self)
            d = self.d;
            i1 = arrayfun(@(i) [i 0 0], 0:d-1, 'UniformOutput', false);
            i2 = arrayfun(@(i) [0 i 0], 0:d-1, 'UniformOutput', false);
            i3 = arrayfun(@(i) [0 0 i], 0:d-1, 'UniformOutput', false);
            D = replab.FiniteGroupDecomposition(self, {i1 i2 i3});
        end
        
        function omega = rootsOfUnity(self)
        % Returns E(d)^0 E(d)^1 ... E(d)^(d-1)
        % where E(d) = exp(2i*pi/d)
            d = self.d;
            switch d
              case 1
                omega = 1;
              case 2
                omega = [1 -1];
              case 4
                omega = [1 1i -1 -1i];
              otherwise
                omega = exp(2i*pi*(0:d-1)/d);
            end
        end
        
        function rep = naturalRep(self)
            omega = self.rootsOfUnity;
            d = self.d;
            W = diag(omega(2)*ones(1, d));
            X = sparse([2:d 1], 1:d, ones(1, d));
            Z = diag(omega);
            if ~replab.Parameters.useSparse
                X = full(X);
            else
                W = sparse(W);
                Z = sparse(Z);
            end
            rep = replab.Rep.lambda(self, 'C', d, true, @(g) W^g(1)*X^g(2)*Z^g(3));
        end

    end

    methods (Access = protected) % Implementation of element enumeration
        
        function el = elementAt(self, index)
            d = self.d;
            ind = index - 1;
            c = double(mod(ind, d));
            ind = (ind - c)/d;
            b = double(mod(ind, d));
            ind = (ind - b)/d;
            a = double(ind);
            el = [a b c];
        end
        
        function index = elementFind(self, x);
            d = self.d;
            a = x(1);
            b = x(2);
            c = x(3);
            index = vpi(a)*d*d + vpi(b)*d + vpi(c + 1);
        end
        
    end

end
