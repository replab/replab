classdef BSGSGroup < replab.FiniteGroup
% Represents a group described by a base and strong generating set
%
% Its construction requires a faithful BSGS action.
    properties
        A; % BSGS action
        G; % Domain derived from the action G = A.G
    end
    
    properties (Access = protected)
        chain_ = [];
        order_ = [];
    end

    methods (Access = protected)
        
        function el = enumeratorAt(self, index)
            indices = self.chain.indicesFromIndex(index);
            el = self.chain.elementFromIndices(indices);
        end
        
        function index = enumeratorFind(self, el)
            [remaining indices] = self.chain.sift(el);
            if isempty(remaining) || ~self.G.isIdentity(remaining)
                ind = [];
            else
                index = self.chain.indexFromIndices(indices);
            end
        end
            
    end
    methods
        
        function E = elements(self)
            E = replab.Enumerator(self.G, self.order, ...
                                  @(ind) self.enumeratorAt(ind), ...
                                  @(el) self.enumeratorFind(el));
        end
        
        function self = BSGSGroup(A, generators, orderOpt)
        % Constructs a BSGS group from a faithful
            assert(isa(A, 'replab.cat.BSGSAction'));
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
            if ~self.knownChain || ~self.chain.areWordsCompleted
                for i = 1:self.nGenerators
                    if self.G.eqv(g, self.generator(i))
                        w = replab.Word.generator(i);
                        return
                    elseif self.G.eqv(g, self.generatorInverse(i))
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
            w = self.chain.factor(g);
        end
        
        function rho = representation(self, d, images, isUnitary)
            T = replab.GeneralLinearGroup(d, false);
            rho = replab.FiniteGroupRep(self, images, isUnitary, T);
        end
        
        function rho = signedPermutationRepresentation(self, d, signedPermutations)
            f = @(g) replab.SignedPermutations(self.domainSize).toMatrix(g);
            images = cellfun(f, signedPermutations, 'UniformOutput', false);
            rho = self.representation(d, images, true);
        end
        
        function rho = permutationRepresentation(self, d, permutations)
            f = @(g) replab.Permutations(self.domainSize).toMatrix(g);
            images = cellfun(f, permutations, 'UniformOutput', false);
            rho = self.representation(d, images, true);
        end
    end

end
