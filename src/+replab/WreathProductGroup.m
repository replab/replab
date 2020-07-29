classdef WreathProductGroup < replab.SemidirectProductGroup
% Wreath product group

    properties (SetAccess = protected)
        n % (integer): Number of copies of the base group
        A % (`.CompactGroup`): Factor of base group
    end

    methods (Static) % WreathProductGroup creation

        function w = make(H, A)
        % Constructs a wreath product group
        %
        % Args:
        %   H (`.PermutationGroup`): Group that permutes the factors
        %   A (`.CompactGroup`): Group whose copies compose the factors
        %
        % Returns:
        %   `.WreathProductGroup`: A specialized instance of `.WreathProductGroup`
            if isa(A, 'replab.FiniteGroup')
                w = replab.prods.WreathProductOfFiniteGroups(H, A);
            else
                w = replab.prods.WreathProductOfCompactGroups(H, A);
            end
        end

    end

    methods

        function self = WreathProductGroup(phi)
            self@replab.SemidirectProductGroup(phi);
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
            p = replab.SymmetricGroup(d*n).compose(ip, basePerm);
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
            for i = 1:ng
                subs{n+1-i} = phiA(base{i});
            end
            p = p(subs{:});
            p = p(:)';
            ip = reshape(1:prod(dims), dims);
            ip = permute(ip, fliplr(n + 1 - h));
            ip = ip(:)';
            p = replab.SymmetricGroup(prod(dims)).compose(ip, p);
        end

    end

    methods % Representations

        function rep = imprimitiveRep(self, Arep)
        % Returns an imprimitive representation of this wreath product
        %
        % It acts on a space of dimension ``self.n * Arep.dimension``, which
        % is a direct sum of copies of ``Arep``. The permutation group acts
        % by permuting the blocks.
        %
        % Args:
        %   Arep (`.Rep`): A representation of `.A`
        %
        % Returns:
        %   `.Rep`: The corresponding imprimitive representation
            rep = replab.wreathproduct.ImprimitiveRep(self, Arep);
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

        function rep = primitiveRep(self, Arep)
        % Returns a primitive representation of this wreath product
        %
        % It acts on a space of dimension ``Arep.dimension^self.n``, which
        % is a tensor product of copies of ``Arep``. The permutation group acts
        % by permuting tensor indices.
        %
        % Args:
        %   Arep (`.Rep`): A representation of `.A`
        %
        % Returns:
        %   `.Rep`: The corresponding primitive representation
            rep = replab.wreathproduct.PrimitiveRep(self, Arep);
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
