classdef NiceFiniteGroup < replab.FiniteSubgroup
    
    properties (SetAccess = protected)
        niceMonomorphism % Injective group homomorphism from this group into a permutation group
                         % Must be valid for elements of the parent group too
    end
    
    properties (Access = protected)
        order_ = [];
        chain_ = [];
    end

    methods
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.FinitelyGeneratedGroup(self);
            if self.knownOrder
                names{1, end+1} = 'order';
                values{1, end+1} = self.order;
            end
        end

        function b = knownOrder(self)
            b = ~isempty(self.order_) || ~isempty(self.chain_);
        end

        function o = order(self)
        % Returns the group order, computing it if necessary
        %
        % Returns:
        %   vpi: Order of this group
            if isempty(self.order_)
                self.order_ = self.chain.order;
            end
            o = self.order_;
        end

        function b = knownChain(self)
        % Tests whether the group BSGS chain has been computed
        %
        % Returns:
        %   logical: True if this group BSGS chain has been computed
            b = ~isempty(self.chain_);
        end
                
        function c = chain(self)
        % Returns the BSGS chain for this group and its nice monomorphism
            if isempty(self.chain_)
                imgId = self.niceMonomorphism(self.identity);
                n = length(imgId);
                S = cellfun(@(x) self.niceMonomorphism(x), self.generators);
                self.chain_ = replab.bsgs1.Chain.makeWithImages(n, S, self, self.generators);
            end
            c = self.chain_;
        end
        
        function rho = rep(self, field, dimension, images)
        % Constructs a finite dimensional representation of this group
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   images (row cell array of matrices): Orthonormal/unitary images of the group generators
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            rho = replab.RepByImages(self, field, dimension, self.niceMonomorphism, images);
        end

        function rho = permutationRep(self, dimension, permutations)
        % Constructs a permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification``
        % to obtain a complex representation.
        %
        % Args:
        %   dimension: Dimension of the representation
        %   permutations (row cell array of permutations): Images of the generators as permutations of size "dimension"
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            S = replab.Permutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, permutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end

        function rho = signedPermutationRep(self, dimension, signedPermutations)
        % Returns a real signed permutation representation of this group
        %
        % The returned representation is real. Use ``rep.complexification``
        % to obtain a complex representation.
        %
        % Args:
        %   dimension: Dximension of the representation
        %   signedPermutations (row cell array of signed permutations): Images of the generators as signed permutations of size "dimension"
        %
        % Returns:
        %   replab.Rep: The constructed group representation
            S = replab.SignedPermutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, signedPermutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end
        
        % CompactGroup methods
        
        function g = sampleUniformly(self)
            g = self.chain.sampleUniformly;
        end
        
        % Finite group methods
        
        function e = elements(self)
        % Returns an enumeration of the group elements
        %
        % Returns:
        %   replab.Enumerator: A space-efficient enumeration of the group elements
            error('Not implemented');
        end
        
        function d = decomposition(self)
        % Returns a decomposition of this group as a product of sets
        %
        % Returns:
        %   replab.FiniteGroupDecomposition: The group decomposition
            error('Not implemented');
        end
               
        function rho = rep(self, field, dimension, images)
        % Constructs a finite dimensional real or complex representation of this group
        %
        %     field: 'R' or 'C' for real or complex
        % dimension: representation dimension
        %    images: 1 x n cell array of matrices providing the images
        %            of the group generators
            rho = replab.RepByImages(self, field, dimension, self.niceMonomorphism, images);
        end

        function rho = permutationRep(self, dimension, permutations)
        % Returns a real permutation representation of this group
        %
        %    dimension: dimension of the representation
        % permutations: row cell array of images of the generators as permutations of size "dimension"
            S = replab.Permutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, permutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end
        
        function rho = signedPermutationRep(self, dimension, signedPermutations)
        % Returns a real signed permutation representation of this group
        %
        %          dimension: dimension of the representation
        % signedPermutations: row cell array of images of the generators as permutations of size "dimension"
            S = replab.SignedPermutations(dimension);
            f = @(g) S.toMatrix(g);
            images = cellfun(f, signedPermutations, 'uniform', 0);
            rho = self.rep('R', dimension, images);
        end

        function rep = leftRegularRep(self)
            o = self.order;
            assert(o < 1e6);
            o = double(o);
            perms = cell(1, self.nGenerators);
            E = self.elements;
            for i = 1:self.nGenerators
                g = self.generator(i);
                img = zeros(1, o);
                for j = 1:o
                    img(j) = double(E.find(self.compose(g, E.at(j))));
                end
                perms{i} = img;
            end
            rep = self.permutationRep(o, perms);
        end

    end

    
    methods
        
        function R = randomBag(self)
            if isequal(self.randomBag_, [])
                self.randomBag_ = replab.RandomBag(self, self.generators);
            end
            R = self.randomBag_;
        end
        
    end

end
