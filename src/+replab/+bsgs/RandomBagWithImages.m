classdef RandomBagWithImages < replab.Str
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
% permutation groups, and keeping track of images of the generated
% elements under some group homomorphism.
%
% When the group homomorphism support is not desired, the `+replab.+bsgs.TrivialGroup`
% trivial group can be used as a placeholder.

    properties (SetAccess = protected)
        n % domainSize
        x0 % Last generated sample
        x % n x r matrix representing the contents of the bag
        J % Image group
        y0 % Image of last generated sample
        y % Images of x
    end

    methods

        function s = headerStr(self)
            s = sprintf('Random bag of %d elements', length(self.x));
        end

        function z = compose(self, x, y)
        % Duplicates self.Permutations(n).compose to avoid reference loops
            z = x(y);
        end

        function xInv = inverse(self, x)
        % Duplicates self.Permutations(n).inverse to avoid reference loops
            xInv = zeros(1, self.n);
            xInv(x) = 1:self.n;
        end

        function [xres yres]  = sample(self)
            r = size(self.x, 2);
            s = randi(r);
            t = randi(r);
            while t == s
                t = randi(r);
            end
            if randi(2) == 2
                if randi(2) == 2 % e = 1
                    self.x(:,s) = self.compose(self.x(:,s)', self.x(:,t)');
                    self.y{s} = self.J.compose(self.y{s}, self.y{t});
                else
                    self.x(:,s) = self.compose(self.x(:,s)', self.inverse(self.x(:,t)'));
                    self.y{s} = self.J.composeWithInverse(self.y{s}, self.y{t});
                end
                self.x0 = self.compose(self.x0, self.x(:,s)');
                self.y0 = self.J.compose(self.y0, self.y{s});
            else
                if randi(2) == 2 % e = 1
                    self.x(:,s) = self.compose(self.x(:,t)', self.x(:,s)');
                    self.y{s} = self.J.compose(self.y{t}, self.y{s});
                else
                    tinv = self.inverse(self.x(:,t)');
                    self.x(:,s) = self.compose(tinv, self.x(:,s)');
                    tinv1 = self.J.inverse(self.y{t});
                    self.y{s} = self.J.compose(tinv1, self.y{s});
                end
                self.x0 = self.compose(self.x(:,s)', self.x0);
                self.y0 = self.J.compose(self.y{s}, self.y0);
            end
            xres = self.x0;
            yres = self.y0;
        end

        function self = RandomBagWithImages(n, generators, r, m, J, images)
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
        %   J (replab.Group, optional): Group structure for images
        %   images (row cell array of elements of ``J``): Images of ``generators``
            self.n = n;
            if nargin < 5
                J = replab.bsgs.TrivialGroup;
                images = arrayfun(@(x) [], 1:length(generators), 'uniform', false);
            end
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
            y = cell(1, r);
            if nGens == 0
                % cater for the special case when generators are empty
                for i = 1:r
                    x(:,i) = 1:n; % identity
                    y{i} = J.identity;
                end
            else
                for i = 1:r
                    ind = mod(i-1, nGens)+1;
                    x(:,i) = generators(:, ind);
                    y{i} = images{ind};
                end
            end
            self.x0 = 1:n; % initially, the identity element
            self.x = x;
            self.J = J;
            self.y0 = J.identity;
            self.y = y;
            for i = 1:m
                self.sample; % perform initial shuffles
            end
        end

    end

end
