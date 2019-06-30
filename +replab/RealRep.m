classdef RealRep < replab.Str
% A finitely generated group real representation on GL_d(R)
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
        centralizerAlgebra_;
        fibers_;
        irreducible_;
    end
    
    methods

        function self = RealRep(group, dimension, images, imagesInv, parent, U, Uinv)
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
            self.M = replab.GroupFun('Mat', @isequal, @() rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
        end

        % Str
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{end+1} = 'images';
            names{end+1} = 'imagesInv';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            for i = 1:self.group.nGenerators
                names{end+1} = sprintf('images(%d)', i);
                values{end+1} = self.images{i};
            end
        end

        function s = shortStr(self, maxColumns)
            if self.isUnitary
                t = 'Symmetric representation';
            else
                t = 'Representation';
            end
            s = sprintf('%s of dimension %d', t, self.dimension);
        end
        
        function lines = longStr(self, maxRows, maxColumns)
            lines = replab.str.longStr(self, maxRows, maxColumns);
            lines{1} = self.shortStr(maxColumns);
        end

        % Own methods
        
        function rr = forget(self)
        % Returns a RealRep that forgot all its special structure
            rr = replab.RealRep(self.group, self.dimension, self.images, self.imagesInv);
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
        
        function b = isMonomial(self)
        % Returns true if this representation is monomial (i.e. is composed of signed permutation matrices)
            b = all(cellfun(@(U) replab.SignedPermutations.isSignedPermutationMatrix(U), self.images));
        end
        
        function b = isTrivial(self)
        % Returns true if this representation is trivial, i.e. all its images are the identity matrix
            b = all(cellfun(@(U) self.M.isIdentity(U), self.images));
        end
        
        function c = centralizerAlgebra(self)
        % Describes the algebra of matrices that commute with this representation
            if isempty(self.centralizerAlgebra_)
                if self.isMonomial
                    self.centralizerAlgebra_ = replab.RealConfigurationAlgebra.fromRealRep(self);
                else
                    self.centralizerAlgebra_ = replab.RealCentralizerAlgebra(self);
                end
            end
            c = self.centralizerAlgebra_;
        end
        
        function sub = slice(self, indices)
        % Returns the subrepresentation given by rho_g(indices, indices), where
        % indices is a subset of the integers {1..self.dimension}
        %
        % The method checks that the given indices correspond to a union of fibers
            bi = unique(self.fibers.blockIndex(indices));
            for b = bi
                assert(all(ismember(self.fibers.blocks{b}, indices)), 'The slice is not a subrepresentation.');
            end
            newImages = cellfun(@(M) M(indices, indices), self.images, 'UniformOutput', false);
            newImagesInv = cellfun(@(M) M(indices, indices), self.imagesInv, 'UniformOutput', false);
            sub = replab.RealRep(self.group, length(indices), newImages, newImagesInv);
        end
        
        function n = nFibers(self)
        % Returns the number of fibers in this representation
            n = self.fibers.nBlocks;
        end
        
        function sub = fiber(self, f)
        % Returns the subrepresentation corresponding to the f-th fiber
            sub = self.slice(self.fibers.blocks{f});
        end
        
        function f = fibers(self)
        % Returns the finest partition of {1..self.dimension} such that the corresponding slices
        % are subrepresentations. This can be used to decompose representations that are already 
        % partly block diagonal.
            if isempty(self.fibers_)
                mask = false(self.dimension, self.dimension);
                for i = 1:length(self.images)
                    im = self.images{i};
                    mask = mask | (abs(im) > replab.Settings.doubleEigTol);
                end
                self.fibers_ = replab.Partition.connectedComponents(mask);
            end
            f = self.fibers_;
        end
        
        function rid = irreducible(self)
        % Returns the decomposition of this representation into irreducible representations
            if isempty(self.irreducible_)
                self.irreducible_ = replab.rep.irreducibleDecomposition(self);
            end
            rid = self.irreducible_;
        end

    end
    
end
