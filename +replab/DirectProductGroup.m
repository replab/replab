classdef DirectProductGroup < replab.FiniteGroup
% Describes an external direct product of finite groups
    
    properties (SetAccess = protected)
        factors; % 1 x nFactors cell array of FiniteGroup instances
    end
    
    methods
        
        function self = DirectProductGroup(factors)
            n = length(factors);
            for i = 1:n
                assert(isa(factors{i}, 'replab.FiniteGroup'), ...
                       'All factors must be finite groups.');
            end
            self.factors = factors;
            identity = cell(1, n);
            generators = {};
            for i = 1:n
                identity{i} = factors{i}.identity;
            end
            for i = 1:n
                factor = factors{i};
                for j = 1:factor.nGenerators
                    generator = identity;
                    generator{i} = factor.generator(j);
                    generators{1, end+1} = generator;
                end
            end
            self.identity = identity;
            self.generators = generators;
        end
        
        function n = nFactors(self)
            n = length(self.factors);
        end
        
        function f = factor(self, i)
            f = self.factors{i};
        end
        
        % Str
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.FiniteGroup(self);
            names{1, end+1} = 'factors';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nFactors
                names{1, end+1} = sprintf('factor(%d)', i);
                values{1, end+1} = self.factor(i);
            end
        end
            
        % Domain
        
        function b = eqv(self, x, y)
            b = true;
            for i = 1:self.nFactors
                if ~self.factor(i).eqv(x{i}, y{i})
                    b = false;
                    return
                end
            end
        end
        
        function g = sample(self)
            g = cell(1, self.nFactors);
            for i = 1:self.nFactors
                g{i} = self.factor(i).sample;
            end
        end
        
        % Semigroup

        function z = compose(self, x, y)
            z = cell(1, self.nFactors);
            for i = 1:self.nFactors
                z{i} = self.factor(i).compose(x{i}, y{i});
            end
        end
        
        % Group
        
        function xInv = inverse(self, x)
            xInv = cell(1, self.nFactors);
            for i = 1:self.nFactors
                xInv{i} = self.factor(i).inverse(x{i});
            end
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, g)
            w = replab.Word.identity;
            shift = 0;
            for i = 1:self.nFactors
                f = self.factor(i);
                wi = f.factorization(g{i});
                w = w * replab.Word.fromIndicesAndExponents(wi.indices + shift, wi.exponents);
                shift = shift + f.nGenerators;
            end
        end
        
        % Finite group
        
        function b = contains(self, g)
            assert(isa(g, 'cell'));
            assert(isvector(g));
            assert(length(g) == self.nFactors);
            b = true;
            for i = 1:self.nFactors
                if ~self.factor(i).contains(g{i})
                    b = false;
                    return
                end
            end
        end
        
        function g = sampleUniformly(self)
            g = cell(1, self.nFactors);
            for i = 1:self.nFactors
                g{i} = self.factor(i).sampleUniformly;
            end
        end
        
        function b = knownOrder(self)
            b = true;
            for i = 1:self.nFactors
                if ~self.factor(i).knownOrder
                    b = false;
                    return
                end
            end
        end
        
        function o = order(self)
            o = vpi(1);
            for i = 1:self.nFactors
                o = o * self.factor(i).order;
            end
        end
        
        function g = atFun(self, ind)
            g = self.identity;
            ind = ind - 1;
            for i = self.nFactors:-1:1
                f = self.factor(i);
                this = mod(ind, f.order);
                ind = (ind - this)/f.order;
                g{i} = f.elements.at(this + 1);
            end
        end
        
        function ind = findFun(self, g)
            ind = vpi(0);
            for i = 1:self.nFactors
                f = self.factor(i);
                ind = ind * f.order;
                ind = ind + f.elements.find(g{i}) - 1;
            end
            ind = ind + 1;
        end
        
        function e = elements(self)
            e = replab.EnumeratorFun(self.order, ...
                                     @(ind) self.atFun(ind), ...
                                     @(g) self.findFun(g));
        end
        
        function gd = decomposition(self)
            T = {};
            for i = 1:self.nFactors
                D = self.factor(i).decomposition;
                Ti = cell(1, length(D));
                for j = 1:length(D)
                    Dj = D{j};
                    Tij = cell(1, length(Dj));
                    for k = 1:length(Dj)
                       Djk = Dj{k};
                       Tijk = self.identity;
                       Tijk{i} = Djk;
                       Tij{k} = Tijk;
                    end
                    Ti{j} = Tij;
                end
                T = horzcat(T, Ti);
            end
            gd = replab.FiniteGroupDecomposition(self, T);
        end
        
    end

end
