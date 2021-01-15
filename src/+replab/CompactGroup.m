classdef CompactGroup < replab.Group
% A group equipped with a Haar measure

    methods % Group construction

        function prd = directProduct(varargin)
        % Returns the direct product of groups
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> G = S3.directProduct(S3, S3);
        %   >>> G.order
        %       216
        %
        % Args:
        %   varargin: Variable number of arguments of type `+replab.CompactGroup`
        %
        % Returns:
        %   `.DirectProductGroup`: The direct product group
        %                                 If all groups are of type `+replab.NiceFiniteGroup`, the return type
        %                                 is `+replab.NiceFiniteGroup` as well.
        %                                 If all groups are of type `+replab.FiniteGroup`, the return type
        %                                 is `+replab.FiniteGroup` as well.
            prd = replab.DirectProductGroup.make(varargin);
        end

        function prd = directPower(self, n)
        % Returns the direct product of this group with itself a number of times
        %
        % Args:
        %   n (integer): Number of copies
        %
        % Returns:
        %   `.CompactGroup`: The direct product self ``x ...(n times)... x self``
        %                    The return type is specialized as in `+replab.CompactGroup.directProduct`.
            factors = arrayfun(@(x) self, 1:n, 'uniform', 0);
            prd = replab.DirectProductGroup.make(factors);
        end

        function sd = semidirectProduct(self, N, phi)
        % Describes an external semidirect product of groups
        %
        % See the construction in https://en.wikipedia.org/wiki/Semidirect_product
        %
        %
        % Let ``H = self`` be a group, ``N`` a group.
        %
        % The semidirect product is defined using a homomorphism
        %
        % `` phi: H -> Aut(N) ``
        %
        % which we write here
        %
        % `` phi: H x N -> N ``, ``n1 = phi(h, n)``
        %
        % Here, we describe this homomorphism by a function handle with parameters of type
        % ``H x N``.
        %
        % We write each semidirect group element {h n}.
        %
        % The type of the return value depends on the most refined type at the intersection
        % of the type of ``self`` and ``N``, with possible types CompactGroup/FiniteGroup/NiceFiniteGroup.
        %
        % Args:
        %   N (`.CompactGroup`): Group acted upon
        %   phi (function_handle): Function describing a homomorphism as described above
        %
        % Returns:
        %   `.CompactGroup`: Semidirect product group
            action = replab.Action.lambda('Semidirect homomorphism', self, N, phi);
            sd = replab.SemidirectProductGroup.make(action);
        end

    end

    methods % Representations

        function rep = trivialRep(self, field, dimension)
        % Returns the trivial representation of this group on a finite dimensional vector space
        %
        % For convenience, either the representation can act on a real or complex vector space,
        % and multiple copies of the 1-dimensional trivial representation can be included, when
        % dimension > 1.
        %
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer, optional): Representation dimension
        %                                  Default value 1
        %
        % Returns:
        %   `.Rep`: An instance of the trivial representation
            if nargin < 2
                dimension = 1;
            end
            rep = replab.rep.TrivialRep(self, field, dimension);
        end

        function rep = directSumRep(self, field, reps)
        % Computes the direct sum of representations on this group
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   reps (cell(1,\*) of `.Rep`): Representation of this group over the given field
        %
        % Returns:
        %   `.Rep`: Direct sum of the representations
            rep = replab.rep.DirectSumRep(self, field, reps);
        end

        function rep = tensorRep(self, field, reps)
        % Computes the tensor product of representations
        %
        % Args:
        %   reps (cell(1,\*) of `.Rep`): Representation of the same group over the same field
        %
        % Returns:
        %   `.Rep`: Tensor product of the representations
            rep = replab.rep.TensorRep(self, field, reps);
        end

    end

    methods (Static) % Group construction

        function group = lambda(header, eqvFun, sampleFun, composeFun, identity, inverseFun)
        % Constructs a compact group from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %   composeFun (function_handle): Handle implementing the `compose` method
        %   identity (element): Identity element of this monoid
        %   inverseFun (function_handle): Handle implementing the `inverse` method
        %
        % Returns:
        %   replab.CompactGroup: The constructed compact group
            group = replab.lambda.CompactGroup(header, eqvFun, sampleFun, composeFun, identity, inverseFun);
        end

    end

end
