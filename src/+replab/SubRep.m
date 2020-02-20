classdef SubRep < replab.Rep
% Describes a subrepresentation of a finite representation
%
% The subrepresentation is described by the basis `H`.

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation
        H_internal % (double(*,*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        F_internal % (double(*,*), may be sparse): Embedding map, dimension ``dParent x dChild``
    end

    methods

        function self = SubRep(parent, H_internal, F_internal)
        % Constructs a subrepresentation of a parent representation
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of which we construct a subrepresentation
        %   H_internal (double(*,*), may be sparse): Subrepresentation basis, dimension ``dParent x dChild``
        %   F_internal (double(*,*), may be sparse): Embedding map, dimension ``dChild x dParent``
            d = size(F_internal, 1);
            dParent = size(F_internal, 2);
            assert(size(H_internal, 1) == dParent);
            assert(size(H_internal, 2) == d);
            assert(parent.dimension == dParent, 'Incorrect basis dimension');
            isUnitary = replab.trileanAnd(parent.isUnitary, isequal(F_internal, H_internal'));
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.isUnitary = isUnitary;
            self.parent = parent;
            self.F_internal = F_internal;
            self.H_internal = H_internal;
        end

        function H = basis(self)
        % Returns the basis of this subrepresentation in the parent representation
        %
        % Always returns a dense matrix.
        %
        % Returns:
        %   double(*,*): Subrepresentation basis given as column vectors
            H = full(self.H_internal);
        end

        %% Str methods

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Rep(self), ...
                {'F_internal' 'H_internal'} ...
                );
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.dimension < 15
                for i = 1:self.dimension
                    names{1, end+1} = sprintf('basis.(:,%d)', i);
                    values{1, end+1} = full(self.H_internal(:,i));
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
% $$$                 X = self.F_internal * self.F_internal';
% $$$                 A = chol(X, 'lower');
% $$$                 U = inv(A) * self.F;
% $$$                 newRep = self.parent.subRepUnitary(U);
% $$$
% $$$         end

        function rho = image_internal(self, g)
            rho = full(self.F_internal*self.parent.image_internal(g)*self.H_internal);
        end

        function rho = inverseImage_internal(self, g)
            rho = full(self.F_internal*self.parent.inverseImage_internal(g)*self.H_internal);
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
            Hs = cellfun(@(sr) sr.H_internal, subReps, 'uniform', 0);
            Fs = cellfun(@(sr) sr.F_internal, subReps, 'uniform', 0);
            newH_internal = horzcat(Hs{:});
            newF_internal = vertcat(Fs{:});
            subRep = parent.subRep(newH_internal, newF_internal);
        end

    end

end
