classdef SemidirectProductGroup < replab.FiniteGroup
% Describes an external semidirect product of finite groups
%
% H is a group, N is a group, and the semidirect product is
% defined using a homomorphism phi: H -> Aut(N)
%
% Here, we describe this homomorphism by an action of the group
% H on the elements of N.
%
% We write each semidirect group element {h n}.
    
    properties (SetAccess = protected)
        H; % group acting
        N; % group acted upon
        phi; % Action of H on N
    end
    
    methods
        
        function self = SemidirectProductGroup(phi)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            self.phi = phi;
            assert(isa(H, 'replab.FiniteGroup'));
            assert(isa(N, 'replab.FiniteGroup'));
            self.H = H;
            self.N = N;
            self.identity = {H.identity N.identity};
            generators = cell(1, H.nGenerators + N.nGenerators);
            for i = 1:length(H.generators)
                generators{i} = {H.generator(i) N.identity};
            end
            for i = 1:length(N.generators)
                generators{H.nGenerators + i} = {H.identity N.generator(i)};
            end
            self.generators = generators;
        end
                    
        % Domain
        
        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end
        
        function g = sample(self)
            g = {self.H.sample self.N.sample};
        end
        
        % Semigroup

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
        
        % Group
        
        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = self.N.inverse(self.phi.leftAction(xh, xn));
            z = {zh zn};
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, g)
            wH = self.H.factorization(g{1});
            wN = self.N.factorization(g{2});
            wNin = replab.Word.fromIndicesAndExponents(wN.indices + self.H.nGenerators, wN.exponents);
            w = wH * wNin;
        end
        
        % Finite group
        
        function b = contains(self, g)
            assert(isa(g, 'cell'));
            assert(isvector(g));
            assert(length(g) == 2);
            b = self.H.contains(g{1}) && self.N.contains(g{2});
        end
        
        function g = sampleUniformly(self)
            g = {self.H.sampleUniformly self.N.sampleUniformly};
        end
        
        function b = knownOrder(self)
            b = self.H.knownOrder && self.N.knownOrder;
        end
        
        function o = order(self)
            o = self.H.order * self.N.order;
        end
        
        function g = atFun(self, ind)
            ind = ind - 1;
            indN = mod(ind, self.N.order);
            indH = (ind - indN)/self.N.order;
            g = {self.H.elements.at(indH + 1) self.N.elements.at(indN + 1)};
        end
        
        function ind = findFun(self, g)
            indH = self.H.elements.find(g{1});
            indN = self.N.elements.find(g{2});
            ind = (indH - 1)*self.N.order + indN;
        end
        
        function e = elements(self)
            e = replab.Enumerator.lambda(self.order, ...
                                         @(ind) self.atFun(ind), ...
                                         @(g) self.findFun(g));
        end
        
        function gd = decomposition(self)
            TH = self.H.decomposition.transversals;
            TN = self.N.decomposition.transversals;
            idN = self.N.identity;
            idH = self.H.identity;
            TH1 = cellfun(@(t) cellfun(@(h) {h idN}, t, 'uniform', 0), TH, 'uniform', 0);
            TN1 = cellfun(@(t) cellfun(@(n) {idH n}, t, 'uniform', 0), TN, 'uniform', 0);
            gd = replab.FiniteGroupDecomposition(self, horzcat(TH1, TN1));
        end
        
    end

end
