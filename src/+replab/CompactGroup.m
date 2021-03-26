classdef CompactGroup < replab.Group
% A group equipped with a Haar measure

    methods % Group properties

        function b = hasReconstruction(self)
        % Returns whether the group has a "reconstruction"
        %
        % Returns:
        %   logical: True if the call to `.reconstruction` would succeed
            b = false;
        end

        function [mu, R] = reconstruction(self)
        % Returns a reconstruction of the group used to speed up group averaging
        %
        % When ``G`` is a compact Lie group, we consider the chain of subgroups ``G >= G0 >= T``, where:
        %
        % - ``G0`` is the connected component of the group around the identity; we also know that the number of left cosets
        %   of ``G0`` in ``G`` is finite; see the group of components, or, for example:
        %   `<https://mathoverflow.net/questions/378160/improved-classification-of-compact-lie-groups>_` .
        %
        %   We write ``R`` a (finite) set of left coset representatives. Thus any element ``g`` of ``G`` can be written ``g = r g0``
        %   where ``r \in R`` and ``g0 \in G0``.
        %
        % - ``G0`` contains a maximal torus ``T``.
        %
        % To average over the action of a representation ``rho`` of ``G`` on ``x``, we write ``x1 = \int_G d\mu(g) \rho(g) x``.
        %
        % However, this is equivalent to ``x1 = 1/|R| \sum_{r \in R} \int_G0 d\mu(g) \rho(g) x``.
        %
        % We finally write ``x1 = 1/{R} \sum_{r \in R} \int_G0 d\mu(g) \rho(g) \int_T d\mu(t) \rho(t) x``; again, the fact that we integrate
        % first over the maximal torus subgroup does not affect the result. However, integration over the maximal torus can often be done
        % very quickly, and the reconstruction speeds up the computation of various equivariant subspaces dramatically.
        %
        % This method returns ``R`` as a `.SetProduct`, and ``T`` as a `.Morphism` from ``T`` to ``G``.
            error('Reconstruction not available');
        end

        function d = maximalTorusDimension(self)
        % Returns, if available, the dimension of the maximal torus contained in the connected component of this group
        %
        % Returns:
        %   integer or ``[]``: Maximal torus dimension if available, or ``[]`` if no reconstruction is available
            if ~self.hasReconstruction
                d = [];
            else
                tm = self.reconstruction;
                d = tm.source.n;
            end
        end

    end

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
        % We write each semidirect group element ``{h n}``.
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

        function rep = commutingRepsRep(self, field, dimension, reps)
        % Constructs a representation from commuting representations
        %
        % Args:
        %   field ({'R', 'C'}): Field
        %   dimension (integer): Dimension of the representation
        %   reps (cell(1,\*) of `.Rep`): Commuting representations with the given field and dimension
        %
        % Returns:
        %   `+replab.Rep`: A representation computed from the product of representations
            rep = replab.rep.CommutingRepsRep(self, field, dimension, reps);
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
