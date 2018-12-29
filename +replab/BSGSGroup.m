classdef BSGSGroup < replab.FiniteGroup
% Represents a group described by a base and strong generating set
%
% Its construction requires a faithful BSGS action.
    properties (SetAccess = protected)
        parent; % Parent group providing composition, inverse, etc...
        action; % faithful action
    end

    properties (Access = protected)
        chain_ = [];
        order_ = [];
    end
    
    methods
        
        function self = BSGSGroup(parent, action, generators, orderOpt)
            assert(isa(action, 'replab.FaithfulAction'));
            self.parent = parent;
            self.action = action;
            self.generators = generators;
            self.identity = parent.identity;
            if nargin > 3 && ~isempty(orderOpt)
                self.order_ = vpi(orderOpt);
            end
        end
        
        % Domain
        
        function b = eqv(self, x, y)
            b = self.parent.eqv(x, y);
        end
        
        function g = sample(self)
            if self.knownChain
                g = self.sampleUniformly;
            else
                g = self.randomBag.sample;
            end
        end
        
        % Semigroup/Monoid/Group
        
        function z = compose(self, x, y)
            z = self.parent.compose(x, y);
        end
        
        function y = inverse(self, x)
            y = self.parent.inverse(x);
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, x)
            if ~self.knownChain || ~self.chain.areWordsCompleted
                for i = 1:self.nGenerators
                    if self.eqv(x, self.generator(i))
                        w = replab.Word.generator(i);
                        return
                    elseif self.eqv(x, self.generatorInverse(i))
                        w = inv(replab.Word.generator(i));
                        return
                    end
                end
                % TODO: tune parameters, inspiration from GAP homomorphisms?
                maxCount = 2000;
                maxLength = max(1, floor(log(maxCount)/log(self.nGenerators)));
                B = length(self.chain.base);
                s = B^2;
                l = B;
                self.chain.wordsQuick(self, maxLength, s, l);
                self.chain.wordsComplete;
            end
            w = self.chain.factor(x);
        end
        
        % FiniteGroup
        
        function b = contains(self, g)
            s = self.chain.sift(g);
            b = ~isequal(s, []) && self.isIdentity(s);
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
        
        function E = elements(self)
            E = replab.EnumeratorFun(self.order, ...
                                     @(ind) self.enumeratorAt(ind), ...
                                     @(el) self.enumeratorFind(el));
        end
        
        function d = decomposition(self)
            d = replab.FiniteGroupDecomposition(self, self.chain.groupDecomposition);
        end
        
        function g = sampleUniformly(self)
            g = self.chain.random;
        end

    end
    
    methods (Access = protected)
        
        function el = enumeratorAt(self, index)
            indices = self.chain.indicesFromIndex(index);
            el = self.chain.elementFromIndices(indices);
        end
        
        function index = enumeratorFind(self, el)
            [remaining indices] = self.chain.sift(el);
            if isempty(remaining) || ~self.parent.isIdentity(remaining)
                ind = [];
            else
                index = self.chain.indexFromIndices(indices);
            end
        end
            
    end
    
    methods
        
        function b = knownChain(self)
            b = ~isempty(self.chain_);
        end
        
        function c = chain(self)
            if isempty(self.chain_)
                self.chain_ = replab.bsgs.Chain.forBSGSGroup(self);
                if isempty(self.order_)
                    self.order_ = self.chain_.order;
                else
                    msg = 'Known group order %s is not equal to computed BSGS chain order %s'; 
                    assert(isequal(self.order_, self.chain_.order), msg, ...
                           num2str(self.order_), num2str(self.chain_.order));
                end
            end
            c = self.chain_;
        end
            
    end

end
