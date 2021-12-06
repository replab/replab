classdef WreathProductGroup < replab.SemidirectProductGroup
% Wreath product group
%
% This group contains:
%
% * a permutation group `.H` acting on $n$ elements,
% * $n$ copies of a compact group `.A` (the base group), collected in the direct product group `.H`,
%
% and those two groups are related through a semidirect product construction (see `.SemidirectProductGroup`), where
% the permutation group acts on `.H` by permuting the copies of `.A`.
%
% Note that compared to Wikipedia `<https://en.wikipedia.org/wiki/Wreath_product>`_ we write the elements of this
% wreath product as a tuple $(h, n)$ where the element $h \in H$ is written first, followed by $n = (a_1, ..., a_n)$
% with $n \in N$ and $a_i \in A$.
%
% The permutation group $H$ acts on $N$ as follows. Let $h \in H$ be a permutation and $n = (a_1, ..., a_n)$, and
% $h(n)$ the result of the action of $h$ on $n$. Then $h(n)_{h(i)} = a_i$, so that $h(n)_i = a_{h^{-1}(i)}$.
%
% According to the literature (see also Wikipedia above), a wreath product defines two canonical actions, the
% imprimitive wreath product action and the primitive wreath product action.
%
% Correspondingly, given a representation of the base group `.A`, we provide methods to construct "imprimitive" and
% "primitive$ representations of this wreath product group using `.imprimitiveRep` and `.primitiveRep`.

    properties (SetAccess = protected)
        n % (integer): Number of copies of the base group
        A % (`.CompactGroup`): Base factor group
    end

    methods (Static) % WreathProductGroup creation

        function w = make(H, A)
        % Constructs a wreath product group
        %
        % Args:
        %   H (`.PermutationGroup`): Complement group that permutes the factors
        %   A (`.CompactGroup`): Base factor group of which the direct product composes the base
        %
        % Returns:
        %   `.WreathProductGroup`: A specialized instance of `.WreathProductGroup`
            if isa(A, 'replab.FiniteGroup')
                type = replab.prods.WreathProductGroupType(H.type, A.type);
                iso = replab.prods.WreathProductNiceIsomorphism(type, H, A);
                w = iso.source;
            else
                w = replab.prods.WreathProductGroup_compact(H, A);
            end
        end

    end

    methods % Permutation actions

        function p = imprimitivePermutation(self, w, phiA)
        % Returns the permutation corresponding to the canonical imprimitive action
        %
        % See https://en.wikipedia.org/wiki/Wreath_product
        %
        % Args:
        %   w (element): Wreath product group element to compute the image of
        %   phi (function_handle, optional): Morphism from elements of A to permutations
        %                                    If omitted default to identity, valid only when
        %                                    A is a permutation group
        %
        % Returns:
        %   permutation: The permutation corresponding to the imprimitive action of ``w``
            if nargin < 3 || isempty(phiA)
                assert(isa(self.A, 'replab.PermutationGroup'));
                phiA = @(x) x;
            end
            n = self.n;
            h = w{1};
            base = w{2};
            im = phiA(base{1});
            d = length(im);
            basePerm = im;
            shift = d;
            for i = 2:n
                im = phiA(base{i}) + shift;
                basePerm = [basePerm im];
                shift = shift + d;
            end
            ip = reshape(1:n*d, [d n]);
            ip = ip(:,h);
            ip = ip(:)';
            S = replab.SymmetricGroup.make(d*n);
            p = S.compose(ip, basePerm);
        end

        function p = primitivePermutation(self, w, phiA)
        % Returns the permutation corresponding to the canonical primitive action
        %
        % See https://en.wikipedia.org/wiki/Wreath_product
        %
        % Args:
        %   w (element): Wreath product group element to compute the image of
        %   phi (function_handle, optional): Morphism from elements of A to permutations
        %                                    If omitted default to identity, valid only when
        %                                    A is a permutation group
        %
        % Returns:
        %   permutation: The permutation corresponding to the primitive action of ``w``
            if nargin < 3 || isempty(phiA)
                assert(isa(self.A, 'replab.PermutationGroup'));
                phiA = @(x) x;
            end
            n = self.n;
            h = w{1};
            base = w{2};
            d = length(phiA(self.A.identity));
            dims = ones(1, n) * d;
            p = reshape(1:prod(dims), dims);
            subs = cell(1, n);
            for i = 1:n
                subs{n+1-i} = phiA(base{i});
            end
            p = p(subs{:});
            p = p(:)';
            ip = reshape(1:prod(dims), dims);
            ip = permute(ip, fliplr(n + 1 - h));
            ip = ip(:)';
            S = replab.SymmetricGroup.make(prod(dims));
            p = S.compose(ip, p);
        end

    end

    methods % Representations

        function rep = imprimitiveRep(self, factorRep)
        % Returns an imprimitive representation of this wreath product
        %
        % It acts on a space of dimension ``self.n * factorRep.dimension``, which
        % is a direct sum of copies of ``factorRep``. The permutation group acts
        % by permuting the blocks.
        %
        % Args:
        %   factorRep (`.Rep`): A representation of the base factor group `.A`
        %
        % Returns:
        %   `.Rep`: The corresponding imprimitive representation
            n = self.n;
            dfr = factorRep.dimension;
            baseRep = self.N.directSumFactorRep(factorRep.field, repmat({factorRep}, 1, n));
            images = cell(1, self.H.nGenerators);
            for i = 1:self.H.nGenerators
                h = self.H.generator(i);
                images{i} = kron(sparse(h, 1:n, ones(1,n), n, n), speye(dfr));
            end
            actingRep = self.H.repByImages(factorRep.field, n*dfr, ...
                                           'preimages', self.H.generators, ...
                                           'images', images);
            rep = self.semidirectProductRep(actingRep, baseRep);
        end

        function rep = imprimitiveRepFun(self, fun)
        % Returns an imprimitive representation of this wreath product
        %
        % See `.imprimitiveRep`.
        %
        % Args:
        %   fun (function_handle): A function that returns a representation of `.A`
        %                          when called on `.A` as in ``Arep = fun(self.A)``
        %
        % Returns:
        %   `.Rep`: The corresponding imprimitive representation
            rep = self.imprimitiveRep(fun(self.A));
        end

        function rep = primitiveRep(self, factorRep)
        % Returns a primitive representation of this wreath product
        %
        % It acts on a space of dimension ``factorRep.dimension^self.n``, which
        % is a tensor product of copies of ``factorRep``. The permutation group acts
        % by permuting tensor indices.
        %
        % Args:
        %   factorRep (`.Rep`): A representation of the base factor group `.A`
        %
        % Returns:
        %   `.Rep`: The corresponding primitive representation
            n = self.n;
            dfr = factorRep.dimension;
            dims = repmat(dfr, 1, n);
            d = prod(dims);
            baseRep = self.N.tensorFactorRep(factorRep.field, repmat({factorRep}, 1, n));
            images = cell(1, self.H.nGenerators);
            for i = 1:self.H.nGenerators
                h = self.H.generator(i);
                I = permute(reshape(1:d, dims), fliplr(n + 1 - h));
                I = I(:)';
                images{i} = sparse(I, 1:d, ones(1, d), d, d);
            end
            actingRep = self.H.repByImages(factorRep.field, dfr^n, ...
                                           'preimages', self.H.generators, ...
                                           'images', images);
            rep = self.semidirectProductRep(actingRep, baseRep);
        end

        function rep = primitiveRepFun(self, fun)
        % Returns an primitive representation of this wreath product
        %
        % See `.primitiveRep`
        %
        % Args:
        %   fun (function_handle): A function that returns a representation of `.A`
        %                          when called on `.A` as in ``Arep = fun(self.A)``
        %
        % Returns:
        %   `.Rep`: The corresponding primitive representation
            rep = self.primitiveRep(fun(self.A));
        end

    end

end
