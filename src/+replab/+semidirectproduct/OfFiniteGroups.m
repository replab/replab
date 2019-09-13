classdef OfFiniteGroups < replab.semidirectproduct.OfCompactGroups & replab.FiniteGroup
% Describes an external semidirect product of finite groups    
    
    methods
        
        function self = OfFiniteGroups(phi)
            self = self@replab.semidirectproduct.OfCompactGroups(phi);
            H = phi.G;
            N = phi.P;
            generators = cell(1, H.nGenerators + N.nGenerators);
            for i = 1:length(H.generators)
                generators{i} = {H.generator(i) N.identity};
            end
            for i = 1:length(N.generators)
                generators{H.nGenerators + i} = {H.identity N.generator(i)};
            end
            self.generators = generators;
        end
                    
        function t = requiredType(self)
            t = 'replab.FiniteGroup';
        end

        %% Domain methods
        
        function g = sample(self)
            g = sample@replab.semidirectproduct.OfCompactGroups(self); % force method selection
        end
        
        %% CompactGroup methods
        
        function g = sampleUniformly(self)
            g = sampleUniformly@replab.semidirectproduct.OfCompactGroups(self); % force method selection
        end
        
        %% FiniteGroup methods
        
        function o = order(self)
            o = self.H.order * self.N.order;
        end
        
        function e = elements(self)
            e = replab.IndexedFamily.lambda(self.order, ...
                                            @(ind) self.atFun(ind), ...
                                            @(g) self.findFun(g));
        end
        
        function gd = decomposition(self)
            TH = self.H.decomposition.T;
            TN = self.N.decomposition.T;
            idN = self.N.identity;
            idH = self.H.identity;
            TH1 = cellfun(@(t) cellfun(@(h) {h idN}, t, 'uniform', 0), TH, 'uniform', 0);
            TN1 = cellfun(@(t) cellfun(@(n) {idH n}, t, 'uniform', 0), TN, 'uniform', 0);
            gd = replab.FiniteGroupDecomposition(self, horzcat(TH1, TN1));
        end

    end

    methods (Access = protected)
        
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
        
    end
    

end
