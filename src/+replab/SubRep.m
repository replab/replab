classdef SubRep < replab.Rep
% Describes a subrepresentation of a finite representation
%
% The subrepresentation is described by the basis `basis`.

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation
        B_internal % (double(\*,\*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        E_internal % (double(\*,\*), may be sparse): Embedding map, dimension ``dParent x dChild``
    end

    methods

        function self = SubRep(parent, B_internal, E_internal)
        % Constructs a subrepresentation of a parent representation
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of which we construct a subrepresentation
        %   B_internal (double(\*,\*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        %   E_internal (double(\*,\*), may be sparse): Embedding map, dimension ``dChild x dParent``
            d = size(E_internal, 1);
            dParent = size(E_internal, 2);
            assert(size(B_internal, 1) == dParent);
            assert(size(B_internal, 2) == d);
            assert(size(E_internal, 1) == d);
            assert(size(E_internal, 2) == dParent);
            assert(parent.dimension == dParent, 'Incorrect basis dimension');
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            if isequal(E_internal, B_internal') && isequal(parent.isUnitary, true)
                self.isUnitary = true;
            end
            self.parent = parent;
            self.E_internal = E_internal;
            self.B_internal = B_internal;
        end

        function verifyInvariance(self)
        % Verifies that the basis of this object defines a subrepresentation up to a given precision
        %
        % From the basis `.basis` which we write $B$ of this subrepresentation, we compute a projector
        % $\pi_R = B B^\dagger$. The user provides a parameter ``epsilon`` ($\epsilon$), and when the verification
        % succeeds, it means that there exists a projector $\pi_K$ on an invariant subspace
        % with $|| \pi_R - \pi_K ||_2 \le \epsilon$.
        %
        % Limitations:
        %
        % * Assumes that the parent representation is unitary, and the basis is
        %   orthonormal.
        % * Assumes that the errors in the basis are much larger than the errors
        %   in the representation images themselves.
        % * Only works for finite groups.
        %
        % Args:
        %   epsilon (double): Maximal error on the subrepresentation basis (see above)
            assert(isa(self.group, 'replab.FiniteGroup'));
            assert(isequal( [isUnitary] )); % ETC
            % by the construction of the BSGS chain, there is a family of sets {T_i},
            % where i = 1,...,depth, such that every group element can be written as
            % g = t_1 t_2 ... t_depth, and t_i \in T_i
            %
            % As S, we take the union of the T_i. Done.
            %
            % Then k = depth, delta = 0, q = 0
            %
            % Use Algorithm IV.1
        end

        function verifyIrreducibilityExact(self)
        % Verifies that this subrepresentation is irreducible
        %
        % * Assumes that the basis is exact
        % * Same assumptions as `.verifyInvariance`
        %
        % Args:
        %   pthr (double): Probability of a false positive (should be ~1e-12 or something)
        %   thetamax (double): Upper bound on theta
            assert(isa(self.group, 'replab.FiniteGroup'));
            % S = the whole group, and we sample from the Haar measure
            % t = depth = 1, choose m such that the theta is << 1
            % quality of the generating set enters in the length of the product
            %
            % Estimate the probability of a false negative given the above
            % (run experimental experiments?)
        end

        function H = basis(self)
        % Returns the basis of this subrepresentation in the parent representation
        %
        % Always returns a dense matrix.
        %
        % Returns:
        %   double(\*,\*): Subrepresentation basis given as column vectors
            H = full(self.B_internal);
        end

        function s = refine(self)
        % Refines the numerical basis of this subrepresentation
        %
        % Returns:
        %   double(\*,\*): Similar subrepresentation with basis precision attempted improvement
            s = replab.rep.refineSubRep(self);
        end

        function [s better] = nice(self)
        % Returns a representation similar to the current subrepresentation, with a nicer basis
        %
        % The "niceness" of the basis is implementation dependent. As of the first implementation
        % of this feature, RepLAB tries to make the basis real, and then with small integer
        % coefficients.
        %
        % The returned subrepresentation is not necessarily unitary.
        %
        % In the case no improvement could be made, the original subrepresentation is returned.
        %
        % Returns:
        %   `+replab.SubRep`: A subrepresentation of ``self.parent``
            s = replab.nice.niceSubRep(self);
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Rep(self), ...
                {'E_internal' 'B_internal'} ...
                );
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.dimension < 15
                for i = 1:self.dimension
                    names{1, end+1} = sprintf('basis.(:,%d)', i);
                    values{1, end+1} = full(self.B_internal(:,i));
                end
            else
                names{1, end+1} = 'basis';
                values{1, end+1} = self.basis;
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.SubRepLaws(self);
        end

        % Rep

% $$$         function [A Ainv] = unitaryChangeOfBasis(self)
% $$$ TODO recover this
% $$$             if isequal(self.parent.isUnitary, true)
% $$$                 X = self.E_internal * self.E_internal';
% $$$                 A = chol(X, 'lower');
% $$$                 U = inv(A) * self.F;
% $$$                 newRep = self.parent.subRepUnitary(U);
% $$$
% $$$         end

        function rho = image_internal(self, g)
            rho = self.E_internal*self.parent.image_internal(g)*self.B_internal;
        end

        function rho = inverseImage_internal(self, g)
            rho = self.E_internal*self.parent.inverseImage_internal(g)*self.B_internal;
        end

    end

    methods (Static)

        function sub = fullSubRep(parent)
        % Creates a full subrepresentation of the given representation, with identity basis
        %
        % Args:
        %   parent (`+replab.Rep`): Representaiton
        %
        % Returns:
        %   `+replab.SubRep`: Subrepresentation identical to ``parent``
            d = parent.dimension;
            sub = parent.subRep(speye(d), speye(d));
            assert(isequal(sub.isUnitary, parent.isUnitary));
            sub.trivialDimension = parent.trivialDimension;
            sub.isIrreducible = parent.isIrreducible;
            sub.frobeniusSchurIndicator = parent.frobeniusSchurIndicator;
            sub.isDivisionAlgebraCanonical = parent.isDivisionAlgebraCanonical;
        end

        function subRep = directSum(parent, subReps)
        % Computes the direct sum of subrepresentations of the same parent representation
        %
        % The subrepresentations must not overlap.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation
        %   subReps (cell(1,\*) of `+replab.SubRep`): A row cell array of subrepresentations of the parent representation
        %
        % Returns:
        %   `+replab.SubRep`: A block-diagonal subrepresentation composed of the given subrepresentations
            Hs = cellfun(@(sr) sr.B_internal, subReps, 'uniform', 0);
            Fs = cellfun(@(sr) sr.E_internal, subReps, 'uniform', 0);
            newB_internal = horzcat(Hs{:});
            newE_internal = vertcat(Fs{:});
            subRep = parent.subRep(newB_internal, newE_internal);
        end

    end

end
