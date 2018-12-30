classdef RealRep < replab.Str
% A finitely generated group real representation on GL_d(R)
    
    properties (SetAccess = protected)
        parent; % either [], or a parent representation of which this is a subrepresentation
                % when parent is not [], the two basis matrices below are defined
        U;    % U has size self.parent.dimension x self.dimension
        Uinv; % Uinv has size self.dimension x self.parent.dimension
              % such that self.image(g) = self.Uinv*self.parent.image(g)*self.U
        group; % Group represented
        dimension; % Representation dimension
        images; % Generator images
        imagesInv; % Generator inverse images
    end
    
    properties (Access = protected)
        M;
        centralizerAlgebra_;
        fibers_;
    end
    
    methods
        
        function self = RealRep(group, dimension, images, imagesInv)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.dimension = dimension;
            self.images = images;
            self.imagesInv = imagesInv;
            d = dimension;
            self.M = replab.GroupFun('Mat', @isequal, @() rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
        end

        function s = str(self)
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

        function b = isUnitary(self)
            b = all(cellfun(@(U) self.M.isIdentity(U*U'), self.images));
        end
        
        function b = isMonomial(self)
            b = all(cellfun(@(U) replab.SignedPermutations.isSignedPermutationMatrix(U), self.images));
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
        
        function I = irreducible(self)
        % Returns the decomposition of this representation into irreducible representations
            I = replab.rep.IsoDec.fromAlgebra(self.centralizerAlgebra);
            I = replab.rep.IrrDec.fromIsoDec(I);
        end

    end
    
end
