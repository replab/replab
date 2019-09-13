classdef OfFiniteGroups < replab.FiniteGroup & replab.directproduct.OfCompactGroups
    
    methods (Access = protected)
        
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

    end
    
    methods
        
        function self = OfFiniteGroups(factors)
            self = self@replab.directproduct.OfCompactGroups(factors);
            generators = {};
            for i = 1:length(factors)
                factor = factors{i};
                for j = 1:factor.nGenerators
                    generator = self.identity;
                    generator{i} = factor.generator(j);
                    generators{1, end+1} = generator;
                end
            end
            self.generators = generators;
        end
        
        %% Str methods

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.directproduct.OfCompactGroups(self), ...
                hiddenFields@replab.FiniteGroup(self) ...
                );
        end
        
        function [names values] = additionalFields(self)
            [names1 values1] = additionalFields@replab.directproduct.OfCompactGroups(self);
            [names2 values2] = additionalFields@replab.FiniteGroup(self);
            names = replab.str.horzcatForce(names1, names2);
            values = replab.str.horzcatForce(values1, values2);
        end
        
        %% Domain methods
        
        function g = sample(self)
            g = sample@replab.directproduct.OfCompactGroups(self); % force method selection
        end
        
        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = sampleUniformly@replab.directproduct.OfCompactGroups(self); % force method selection
        end

        %% FiniteGroup methods
        
        function o = order(self)
            o = vpi(1);
            for i = 1:self.nFactors
                o = o * self.factor(i).order;
            end
        end
                
        function e = elements(self)
            e = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.atFun(ind), ...
                                            @(g) self.findFun(g));
        end
        
        function gd = decomposition(self)
            T = {};
            for i = 1:self.nFactors
                D = self.factor(i).decomposition.T;
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
