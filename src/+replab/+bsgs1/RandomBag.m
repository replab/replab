% Generator of random elements from generators of a permutation group.
% 
% A random bag is a set of random group elements that always generates
% the group; random elements are provided by multiplying elements of the
% bag and returning one element of the product which is removed from the bag.
% 
% Straight-forward implementation of PRINITIALIZE and PRRANDOM of 
% section 3.2.2, pp. 70-71 of Holt 2005 (Handbook of Computational Group Theory)
classdef RandomBag < replab.Str
    
    properties (SetAccess = protected)
        G % Group definition
        x0 % Last generated sample
        x % n x r matrix representing the contents of the bag
        withImages % (logical) Whether we store images as well
        J % Image group
        y0 % Image of last generated sample
        y % Images of x
    end

    methods
        
        function s = headerStr(self)
            s = sprintf('Random bag of %d elements', length(self.x));
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
                    self.x(:,s) = self.G.compose(self.x(:,s)', self.x(:,t)');
                    if self.withImages
                        self.y{s} = self.J.compose(self.y{s}, self.y{t});
                    end
                else
                    self.x(:,s) = self.G.composeWithInverse(self.x(:,s)', self.x(:,t)');
                    if self.withImages
                        self.y{s} = self.J.composeWithInverse(self.y{s}, self.y{t});
                    end
                end
                self.x0 = self.G.compose(self.x0, self.x(:,s)');
                if self.withImages
                    self.y0 = self.J.compose(self.y0, self.y{s});
                end
            else
                if randi(2) == 2 % e = 1
                    self.x(:,s) = self.G.compose(self.x(:,t)', self.x(:,s)');
                    if self.withImages
                        self.y{s} = self.J.compose(self.y{t}, self.y{s});
                    end
                else
                    tinv = self.G.inverse(self.x(:,t)');
                    self.x(:,s) = self.G.compose(tinv, self.x(:,s)');
                    if self.withImages
                        tinv1 = self.J.inverse(self.y{t});
                        self.y{s} = self.J.compose(tinv1, self.y{s});
                    end
                end
                self.x0 = self.G.compose(self.x(:,s)', self.x0);
                if self.withImages
                    self.y0 = self.J.compose(self.y{s}, self.y0);
                end
            end
            xres = self.x0;
            if self.withImages
                yres = self.y0;
            end
        end
        
        function self = RandomBag(n, generators, r, m, J, images)
        % Constructs a random bag from the given permutations
        %
        % Args:
        %   n: Domain size
        %   generators (integer matrix): Group generators, given as a n x nGens matrix
        %   r (integer, optional): Number of elements in the bag
        %                          Must be >= nGens and >= 10
        %                          Default value is max(nGens, 10)
        %   m (integer, optional): Number of shuffles done during initialization 
        %                          Default value is 50
        %  
            G = replab.Permutations(n);
            if nargin < 5
                self.withImages = false;
            else
                self.withImages = true;
                self.J = J;
            end
            if nargin < 5 || isequal(r, [])
                r = -1;
            end
            if nargin < 4 || isequal(m, [])
                m = 50;
            end
            nGens = size(generators, 2); % number of generators
            if r < nGens || r < 10
                r = max(nGens, 10);
            end
            x = zeros(n, r);
            if self.withImages
                y = cell(1, r);
            end
            if nGens == 0
                % cater for the special case when generators are empty
                for i = 1:r
                    x(:,i) = 1:n;
                    if self.withImages
                        y{i} = J.identity;
                    end
                end
            else
                for i = 1:r
                    ind = mod(i-1, nGens)+1;
                    x(:,i) = generators(:, ind);
                    if self.withImages
                        y{i} = images{ind};
                    end
                end
            end
            self.G = G;
            self.x0 = G.identity; % initially, the identity element
            self.x = x;
            if self.withImages
                self.y0 = J.identity;
                self.y = y;
            end
            for i = 1:m
                self.sample; % perform initial shuffles
            end
        end
   
    end
    
end
