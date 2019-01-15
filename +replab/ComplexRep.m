classdef ComplexRep < replab.Str
% A finitely generated group real representation on GL_d(C)
%
% It optionally keeps track of itself being a subrepresentation of a larger parent
% representation
    properties (SetAccess = protected)
        parent; % Either [], or a parent representation of which this is a subrepresentation
                % when parent is not [], the two basis matrices below are defined
        
        U;    % U has size self.parent.dimension x self.dimension
        Uinv; % Uinv has size self.dimension x self.parent.dimension
              % such that self.image(g) = self.Uinv*self.parent.image(g)*self.U
              %
              % Note that we are not using parent in a recursive way: a subrepresentation
              % of a subrepresentation will have its parent refering to the most general
              % representation present. Also said: if self.parent is not [], we still have
              % self.parent.parent = []
        
        group; % Group represented
        dimension; % Representation dimension
        images; % Generator images
        imagesInv; % Generator inverse images
    end
    
    properties (Access = protected)
        M;
    end
    
    methods
        
        function self = ComplexRep(group, dimension, images, imagesInv, parent, U, Uinv)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            if nargin < 5
                parent = [];
                U = [];
                Uinv = [];
            end
            self.group = group;
            self.dimension = dimension;
            self.images = images;
            self.imagesInv = imagesInv;
            self.parent = parent;
            self.U = U;
            self.Uinv = Uinv;
            d = dimension;
            self.M = replab.GroupFun('Mat', @isequal, @() rand(d, d)+1i*rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
        end

        function rr = forget(self)
        % Returns a ComplexRep that forgot all its special structure
            rr = replab.ComplexRep(self.group, self.dimension, self.images, self.imagesInv);
        end

        function s = str(self)
        % Nice string representation
            if self.isUnitary
                t = 'Unitary representation';
            else
                t = 'Representation';
            end
            s = sprintf('%s of dimension %d with generator images', t, self.dimension);
            for i = 1:length(self.images)
                gen = char('a' + i - 1);
                s = [s char(10) '- ' gen ':' char(10)];
                img = replab.prependLines(replab.strOf(self.images{i}), '    ');
                s = [s img char(10)];
            end
        end
                
        function rho = image(self, g)
        % Computes the image of a group element g in this representation
            word = self.group.factorization(g);
            rho = self.M.identity;
            for i = 1:length(word.indices)
                ind = word.indices(i);
                e = word.exponents(i);
                if e > 0
                    ge = self.M.composeN(self.images{ind}, e);
                else
                    ge = self.M.composeN(self.imagesInv{ind}, -e);
                end
                rho = self.M.compose(rho, ge);
            end
        end
        
        function rho = sample(self)
        % Returns the representation image of a random group element
            rho = self.image(self.group.sample);
        end

        function b = isUnitary(self)
        % Returns true if this representation is unitary
            b = all(cellfun(@(U) self.M.isIdentity(U*U'), self.images));
        end
        
    end
    
end
