classdef RandomBag < replab.Str
% Generator of random elements from generators of a permutation group.
%
% A random bag is a set of random group elements that always generates
% the group; random elements are provided by multiplying elements of the
% bag and returning one element of the product which is removed from the bag.
%
% Straight-forward implementation of PRINITIALIZE and PRRANDOM of
% section 3.2.2, pp. 70-71 of Holt 2005 (Handbook of Computational Group Theory)
%
% This implementation differs from `+replab.RandomBag` by specializing for
% permutation groups.

    properties (SetAccess = protected)
        n % (integer): Size of the domain
        x0 % (integer(n, 1)): Last generated sample
        x % (integer(n, \*)): Matrix representing the contents of the bag
    end

    methods

        function s = headerStr(self)
            s = sprintf('Random bag of %d elements', length(self.x));
        end

        function z = compose(self, x, y)
        % Duplicates replab.PermutationGroup.compose to avoid reference loops
            z = x(y);
        end

        function xInv = inverse(self, x)
        % Duplicates replab.PermutationGroup.inverse to avoid reference loops
            xInv = zeros(1, self.n);
            xInv(x) = 1:self.n;
        end

        function xres  = sample(self)
            r = size(self.x, 2);
            s = randi(r);
            t = randi(r);
            while t == s
                t = randi(r);
            end
            xs = self.x(:, s);
            xt = self.x(:, t);
            type = randi(4);
            if type <= 2
                if type == 1
                    xs = xs(xt); % xs = compose(xs, xt)
                else
                    xs(xt) = xs; % xs = compose(xs, inverse(xt))
                end
                self.x0 = self.x0(xs); % x0 = compose(x0, xs)
            else
                if type == 3
                    xs = xt(xs); % xs = compose(xt, xs)
                else
                    % computes xs = compose(inverse(xt), xs);
                    xt(xt) = xt; % invert xt
                    xs = xt(xs); % compose
                end
                self.x0 = xs(self.x0); % x0 = compose(xs, x0)
            end
            self.x(:, s) = xs; % store new xs
            xres = self.x0';
        end

        function self = RandomBag(n, generators, r, m)
        % Constructs a random bag from the given permutations
        %
        % If a image group ``J`` and generator images are not given, the trivial group is used.
        %
        % Args:
        %   n: Domain size
        %   generators (integer matrix): Group generators, given as a n x nGens matrix
        %   r (integer, optional): Number of elements in the bag
        %                          Must be >= nGens and >= 10
        %                          Default value is max(nGens, 10)
        %   m (integer, optional): Number of shuffles done during initialization
        %                          Default value is 50
            self.n = n;
            if nargin < 4 || isempty(m)
                m = 50;
            end
            if nargin < 3 || isempty(r)
                r = -1;
            end
            nGens = size(generators, 2); % number of generators
            if r < nGens || r < 10
                r = max(nGens, 10);
            end
            x = zeros(n, r);
            if nGens == 0
                % cater for the special case when generators are empty
                for i = 1:r
                    x(:,i) = 1:n; % identity
                end
            else
                for i = 1:r
                    ind = mod(i-1, nGens)+1;
                    x(:,i) = generators(:, ind);
                end
            end
            self.x0 = 1:n; % initially, the identity element
            self.x = x;
            for i = 1:m
                self.sample; % perform initial shuffles
            end
        end

    end

end
