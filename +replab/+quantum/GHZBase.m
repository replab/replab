classdef GHZBase < replab.FiniteGroup
% A discrete approximation of the GHZ abelian invariant group

    properties
        rootOrder;
        nParties;
        nLevels;
    end
    
    methods
        
        function self = GHZBase(nParties, nLevels, rootOrder)
            self.nParties = nParties;
            self.nLevels = nLevels;
            self.rootOrder = rootOrder;
            self.identity = zeros(nParties, nLevels);
            generators = cell(1, (nParties-1)*(nLevels-1));
            ind = 1;
            for i = 1:nParties-1
                for j = 2:nLevels
                    g = zeros(nParties, nLevels);
                    g(i,j) = 1;
                    g(nParties,j) = rootOrder - 1;
                    generators{ind} = g;
                    ind = ind + 1;
                end
            end
            self.generators = generators;
        end
        
        function g1 = canonical(self, g)
            ro = self.rootOrder;
            g1 = g;
            for i = 1:self.nParties
                s = ro - g(i,1);
                g1(i,:) = mod(g(i,:) + s, ro);
            end
        end
        
        function p = partial(self, g)
            p = g(1:end-1, 2:end);
        end
        
        function g = fromPartial(self, p)
            nP = self.nParties;
            nL = self.nLevels;
            ro = self.rootOrder;
            last = mod(ro - mod(sum(p, 1), ro), ro);
            g = [zeros(nP-1, 1) p
                 0 last];
        end
        
        function rho = toMatrix(self, g)
            rho = 1;
            for i = 1:self.nParties
                D = diag(exp(1i*2*pi*g(i,:)/self.rootOrder));
                rho = kron(rho, D);
            end
        end
        
        function rep = naturalRep(self)
            d = self.nLevels^self.nParties;
            rep = replab.RepFun(self, 'C', d, @(g) self.toMatrix(g));
        end
        
        function g1 = permuteParties(self, p, g)
            g1(p, :) = g;
        end
        
        function g1 = permuteLevels(self, p, g)
            g1(:, p) = g;
            g1 = self.canonical(g1);
        end
        
        % Domain
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function g = sample(self)
            g = randi(self.rootOrder, self.nParties-1, self.nLevels-1);
            g = self.fromPartial(g - 1);
        end
        
        % Semigroup
        
        function z = compose(self, x, y)
            z = mod(x + y, self.rootOrder);
        end
        
        % Group
        
        function z = inverse(self, x)
            z = mod(self.rootOrder - x, self.rootOrder);
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, g)
            g = g(1:end-1, 2:end)';
            g = g(:);
            n = length(g);
            indices = [];
            exponents = [];
            for i = 1:length(g)
                if g(i) > 0
                    indices = [indices i];
                    exponents = [exponents g(i)];
                end
            end
            w = replab.Word.fromIndicesAndExponents(indices, exponents);
        end
        
        % FiniteGroup
        
        function b = contains(self, g)
            assert(isa(g, 'double'));
            assert(ismatrix(g));
            assert(size(g, 1) == self.nParties);
            assert(size(g, 2) == self.nLevels);
            assert(all(all(g >= 0)));
            assert(all(all(g < self.rootOrder)));
            b = true;
        end
        
        function g = sampleUniformly(self)
            g = self.sample;
        end
           
        function b = knownOrder(self)
            b = true;
        end
        
        function o = order(self)
            n = (self.nParties - 1) * (self.nLevels - 1);
            o = vpi(self.rootOrder)^n;
        end
        
        function g = atFun(self, ind)
            p = zeros(1, (self.nParties-1)*(self.nLevels-1));
            ind = ind - 1;
            for i = length(p):-1:1
                this = double(mod(ind, self.rootOrder));
                ind = (ind - this)/self.rootOrder;
                p(i) = this;
            end
            p = reshape(p, [self.nLevels-1 self.nParties-1])';
            g = self.fromPartial(p);
        end
        
        function ind = findFun(self, g)
            p = self.partial(g);
            p = p';
            p = p(:);
            ind = vpi(0);
            for i = 1:length(p)
                ind = ind * self.rootOrder;
                ind = ind + p(i);
            end
            ind = ind + 1;
        end
        
        function e = elements(self)
            e = replab.EnumeratorFun(self.order, ...
                                     @(ind) self.atFun(ind), ...
                                     @(g) self.findFun(g));
        end
        
        function gd = decomposition(self)
            F = factor(self.rootOrder);
            T = {};
            for i = 1:self.nGenerators
                g = self.generator(i);
                p = 1;
                for j = 1:length(F)
                    t = {self.identity};
                    for k = 1:F(j)-1
                        t{k+1} = g*k*p;
                    end
                    T{1,end+1} = t;
                    p = p * F(j);
                end
            end
            gd = replab.FiniteGroupDecomposition(self, T);
        end
        
    end
    
end
