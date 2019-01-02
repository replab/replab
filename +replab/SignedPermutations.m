classdef SignedPermutations < replab.SignedPermutationGroup
% Describes the signed permutation group over
% {-n...-1, 1...n} where n = domainSize
    
    methods
        
        function self = SignedPermutations(domainSize)
            self = self@replab.SignedPermutationGroup(domainSize);
            switch self.domainSize
              case 0
                self.generators = cell(1, 0); % TODO: verify
              case 1
                self.generators = {[-1]};
              case 2
                self.generators = {[2 1] [-1 2]};
              otherwise
                shift = [2:self.domainSize 1];
                trans = [2 1 3:self.domainSize];
                flip = [-1 2:self.domainSize];
                self.generators = {shift trans flip};
            end
        end
        
        % Str
                
        function s = str(self)
            s = sprintf('Signed permutations acting on {-%d..-1 1..%d}', self.domainSize, self.domainSize);
        end
        
        % Domain
        
        function s = sample(self)
            n = self.domainSize;
            s = randperm(n) .* (randi([0 1], 1, n)*2-1);
        end
        
        % FinitelyGeneratedGroup
        
        function w = factorization(self, x)
            if self.isIdentity(x)
                w = replab.Word.identity;
                return
            elseif self.domainSize == 1
                % not identity, so it is the flip
                w = replab.Word.generator(1);
                return
            elseif self.domainSize == 2
                t = replab.Word.generator(1);
                f1 = replab.Word.generator(2);
                f2 = t * f1 * t;
                w = replab.Word.identity;
                if x(1) < 0
                    x(1) = abs(x(1));
                    w = f1 * w;
                end
                if x(2) < 0
                    x(2) = abs(x(2));
                    w = f2 * w;
                end
                if isequal(x, [2 1])
                    w = t * w;
                end
                return
            end
            n = length(x);
            w = replab.Word.identity;
            % Flip negative images
            for i = 1:n
                if x(i) < 0
                    if i == 1
                        shift = replab.Word.identity;
                    else
                        shift = replab.Word.fromIndicesAndExponents(1, i - 1);
                    end
                    x(i) = abs(x(i));
                    w = shift * replab.Word.generator(3) * inv(shift) * w;
                end
            end
            % Bubble sort to order the image
            moved = true;
            while moved
                moved = false;
                for i = 1:n-1
                    if x(i) > x(i+1)
                        t = x(i+1);
                        x(i+1) = x(i);
                        x(i) = t;
                        moved = true;
                        if i == 1
                            shift = replab.Word.identity;
                        else
                            shift = replab.Word.fromIndicesAndExponents(1, i - 1);
                        end
                        w = shift * replab.Word.generator(2) * inv(shift) * w;
                    end
                end
            end
        end
        
        % Contains
        
        function b = contains(self, g)
            b = length(g) == self.domainSize;
        end
        
        function b = knownOrder(self)
            b = true;
        end
        
        function o = order(self)
            o = factorial(vpi(self.domainSize)) * vpi(2)^self.domainSize;
        end        
        
        function E = elements(self)
            E = replab.EnumeratorFun(self.order, ...
                                     @(ind) self.enumeratorAt(ind), ...
                                     @(el) self.enumeratorFind(el));
        end

        function d = decomposition(self)
            G = self.subgroup(self.generators, self.order);
            d = G.decomposition;
        end
        
        function G = subgroup(self, generators, orderOpt)
            if nargin < 3
                orderOpt = [];
            end
            G = replab.SignedPermutationSubgroup(self, generators, orderOpt);
        end
        
    end
    
    methods (Access = protected)
        
        function ind = enumeratorFind(self, g)
            n = self.domainSize;
            ind0 = vpi(0);
            els = [-n:-1 1:n];
            for i = 1:n
                ind0 = ind0 * 2*(n - i + 1);
                ind0 = ind0 + (find(els == g(i)) - 1);
                els = setdiff(els, [g(i) -g(i)]);
            end
            ind = ind0 + 1;
        end
        
        function g = enumeratorAt(self, ind)
            n = self.domainSize;
            ind0 = ind - 1; % make it 0-based
            inds = zeros(1, n);
            for i = 1:n
                r = mod(ind0, 2*i);
                ind0 = (ind0 - r)/(2*i);
                inds(i) = double(r + 1);
            end
            inds = fliplr(inds);
            els = [-n:-1 1:n];
            g = zeros(1, n);
            for i = 1:n
                e = els(inds(i));
                g(i) = e;
                els = setdiff(els, [e -e]);
            end
        end
        
    end
    
end
