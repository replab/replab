classdef SignedSymmetricGroup < replab.SignedPermutationGroup
% Describes the signed permutation group over {-n...-1, 1...n} where n = domainSize

    methods

        function self = SignedSymmetricGroup(domainSize, genericIsomorphism)
        % Constructs the group of all signed permutations over a given domain size
        %
        % Args:
        %   domainSize (integer): Domain size, must be > 0
            generators = {[-1 2:domainSize]};
            if domainSize > 1
                generators{1,end+1} = [2:domainSize 1];
            end
            if domainSize > 2
                generators{1,end+1} = [2 1 3:domainSize];
            end
            type = 'self';
            self@replab.SignedPermutationGroup(domainSize, generators, 'type', type);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = sprintf('Signed permutations acting on {-%d..-1 1..%d}', self.domainSize, self.domainSize);
        end

        % Domain

        function s = sample(self)
            n = self.domainSize;
            s = randperm(n) .* (randi([0 1], 1, n)*2-1);
        end

        % FiniteGroup

        function b = contains(self, g)
            assert(length(g) == self.domainSize, 'Signed permutation in wrong domain');
            assert(all(g ~= 0) && all(abs(g) <= self.domainSize), 'Signed permutation has out of range coefficients');
            b = true;
        end

    end

    methods (Access = protected)

        function o = computeOrder(self)
            ds = self.domainSize;
            o = replab.util.factorial(ds)*replab.util.multiplyIntegers(ones(1, ds)*2);
            % faster than o = factorial(vpi(domainSize))*vpi(2)^domainSize;
        end

% $$$         function E = computeElementsSequence(self)
% $$$             E = replab.Sequence.lambda(self.order, ...
% $$$                                             @(ind) self.enumeratorAt(ind), ...
% $$$                                             @(el) self.enumeratorFind(el));
% $$$         end

    end

    methods

% $$$         function ind = enumeratorFind(self, g)
% $$$             n = self.domainSize;
% $$$             ind0 = vpi(0);
% $$$             els = [-n:-1 1:n];
% $$$             for i = 1:n
% $$$                 ind0 = ind0 * 2*(n - i + 1);
% $$$                 ind0 = ind0 + (find(els == g(i)) - 1);
% $$$                 els = setdiff(els, [g(i) -g(i)]);
% $$$             end
% $$$             ind = ind0 + 1;
% $$$         end
% $$$
% $$$         function g = enumeratorAt(self, ind)
% $$$             n = self.domainSize;
% $$$             ind0 = ind - 1; % make it 0-based
% $$$             inds = zeros(1, n);
% $$$             for i = 1:n
% $$$                 r = mod(ind0, 2*i);
% $$$                 ind0 = (ind0 - r)/(2*i);
% $$$                 inds(i) = double(r + 1);
% $$$             end
% $$$             inds = fliplr(inds);
% $$$             els = [-n:-1 1:n];
% $$$             g = zeros(1, n);
% $$$             for i = 1:n
% $$$                 e = els(inds(i));
% $$$                 g(i) = e;
% $$$                 els = setdiff(els, [e -e]);
% $$$             end
% $$$         end

    end

end
