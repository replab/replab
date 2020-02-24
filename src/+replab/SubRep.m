classdef SubRep < replab.Rep
% Describes a subrepresentation of a finite representation
%
% The subrepresentation is described by the basis `H`.

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation
        B_internal % (double(*,*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        E_internal % (double(*,*), may be sparse): Embedding map, dimension ``dParent x dChild``
    end

    methods

        function self = SubRep(parent, B_internal, E_internal)
        % Constructs a subrepresentation of a parent representation
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of which we construct a subrepresentation
        %   B_internal (double(*,*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        %   E_internal (double(*,*), may be sparse): Embedding map, dimension ``dChild x dParent``
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
            isUnitary = replab.trileanAnd(parent.isUnitary, isequal(E_internal, B_internal'));
            self.isUnitary = isUnitary;
            self.parent = parent;
            self.E_internal = E_internal;
            self.B_internal = B_internal;
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
            s = replab.rep.refine(self);
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

        %% Str methods

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

        %% Rep methods

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
            rho = full(self.E_internal*self.parent.image_internal(g)*self.B_internal);
        end

        function rho = inverseImage_internal(self, g)
            rho = full(self.E_internal*self.parent.inverseImage_internal(g)*self.B_internal);
        end

    end

    methods (Static)

        function subRep = directSum(parent, subReps)
        % Computes the direct sum of subrepresentations of the same parent representation
        %
        % The subrepresentations must not overlap.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation
        %   subReps (cell(1,*) of `+replab.SubRep`): A row cell array of subrepresentations of the parent representation
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
