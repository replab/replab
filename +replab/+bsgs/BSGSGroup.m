classdef BSGSGroup < replab.FiniteGroup
    
    properties
        G; % Domain
        A; % BSGS action
    end
    
    properties (Access = protected)
        chain_ = [];
        order_ = [];
    end

    methods
        
        function self = BSGSGroup(A, generators, orderOpt)
            self.A = A;
            self.G = A.G;
            self.generators = generators;
            if ~isempty(orderOpt)
                self.order_ = vpi(orderOpt);
            end
        end

        function b = contains(self, g)
            s = self.chain.sift(g);
            b = ~isequal(s, []) && self.G.isIdentity(s);
        end

        function b = knownOrder(self)
            b = ~isempty(self.order_);
        end

        function o = order(self)
            if isempty(self.order_)
                self.order_ = self.chain.order;
            end
            o = self.order_;
        end

        function b = knownChain(self)
            b = ~isempty(self.chain_);
        end
        
        function c = chain(self)
            if isempty(self.chain_)
                self.chain_ = replab.bsgs.Chain.forBSGSGroup(self);
            end
            c = self.chain_;
        end
        
        function g = sample(self)
            if self.knownChain
                g = self.sampleUniformly;
            else
                g = self.randomBag.sample;
            end
        end
        
        function g = sampleUniformly(self)
            g = self.chain.random;
        end
        
        function w = factorization(self, g)
            if ~self.chain.areWordsCompleted
                maxCount = 2000;
                maxLength = max(1, floor(log(maxCount)/log(self.nGenerators)));
                B = length(self.chain.base);
                s = B^2;
                l = B;
                self.chain.wordsQuick(self, maxLength, s, l);
                self.chain.wordsComplete;
            end
            w = self.chain.factor(g);
        end
        
    end

end
