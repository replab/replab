classdef SignedPermutationSubgroup < replab.SignedPermutationGroup & replab.FiniteGroup
% Represents a group of signed permutations
    
    properties (Access = protected)
        nice_ = [];
    end
    
    methods

        function self = SignedPermutationSubgroup(parent, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            self@replab.SignedPermutationGroup(parent.domainSize);
            self.generators = generators;
            if nargin > 2 && ~isempty(orderOpt)
                self.order_ = vpi(orderOpt);
            end
        end
        
        function g = sample(self)
            if self.nice.knownChain
                g = self.fromPermutation(self.nice.sampleUniformly);
            else
                g = self.randomBag.sample;
            end
        end
        
        function G = nice(self)
        % "Nice" permutation group isomorphic to this subgroup
            if isempty(self.nice_)
                gens1 = cellfun(@(x) self.toPermutation(x), self.generators, 'UniformOutput', false);
                self.nice_ = replab.Permutations(self.domainSize * 2).subgroup(gens1);
            end
            G = self.nice_;                
        end
        
        function w = factorization(self, x)
            w = self.nice.factorization(self.toPermutation(x));
        end
        
        function b = contains(self, x)
            b = self.nice.contains(self.toPermutation(x));
        end
        
        function b = knownOrder(self)
            b = self.nice.knownOrder;
        end
        
        function o = order(self)
            o = self.nice.order;
        end
        
        function E = elements(self)
            ea = @(ind) self.fromPermutation(self.nice.elements.at(ind));
            ef = @(x) self.nice.elements.find(self.toPermutation(x));
            E = replab.EnumeratorFun(self.order, ea, ef);
        end
        
        function d = decomposition(self)
            d1 = self.nice.decomposition;
            T1 = d1.transversals;
            T = cellfun(@(t) cellfun(@(x) self.fromPermutation(x), t, 'UniformOutput', false), T1, 'UniformOutput', false);
            d = replab.FiniteGroupDecomposition(self, T);
        end
            
        function g = sampleUniformly(self)
            g = self.fromPermutation(self.nice.sampleUniformly);
        end
        
    end
        
end
