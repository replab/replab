classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
%
% The basis of this subrepresentation in the parent representation is formed by the rows of the
% matrix ``U = D0 * U0``, and we have ``image(g) = U * parent.image(g) * U'``,
% where U0 is an orthogonal basis, not necessarily normalized, and D0 is a correction factor diagonal matrix.

    properties (SetAccess = protected)
        parent % replab.Rep: Parent representation
        U0 % double matrix, can be sparse: Orthogonal basis of dimension dChild x dParent
        D0 % double sparse matrix: Sparse diagonal matrix of correction factors of dimension dChild x dChild
        hasCorrection % logical: True when D0 differs from the identity matrix
        extra % struct: Contains key/value pairs concerning the irreducible decomposition process
              %         'hasTrivialSubspace' whether this contains a trivial subrepresentation
              %         'reducedBlocks' whether the block diagonal reduction heuristic has been applied
              %         'isIrreducible' whether the representation is irreducible
              %         'divisionAlgebra', if known and self.field == 'R', is the division algebra type
              %                            which can be 'R', 'C', or 'H'
              %         'isDivisionAlgebraCanonical', if divisionAlgebra is 'C' or 'H', whether
              %                                       the representation is expressed in the RepLAB
              %                                       canonical basis for that algebra
    end
    
    methods
        
        function self = SubRep(parent, U0, extra)
        % Constructs a subrepresentation of a parent representation
        %
        % Args:
        %   parent (replab.Rep): Parent representation of which we construct a subrepresentation
        %   U0 (double matrix, can be sparse): Basis matrix of dimension dThisChild x dParent
            if nargin < 3
                extra = struct;
            end
            d = size(U0, 1);
            dParent = size(U0, 2);
            assert(parent.dimension == dParent);
            assert(isequal(parent.isUnitary, true), 'We support only unitary representations');
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.isUnitary = true;
            self.extra = extra;
            self.U0 = U0;
            D0 = speye(d);
            hasCorrection = false;
            for i = 1:d
                nrm = self.U0(i,:) * self.U0(i,:)';
                if abs(nrm - 1) > replab.Settings.doubleEigTol
                    D0(i,i) = 1/sqrt(nrm);
                    hasCorrection = true;
                end
            end
            self.D0 = D0;
            self.hasCorrection = hasCorrection;
        end
        
        function b = isExtraTrue(self, name)
        % Returns whether the value of the given extra is known AND true
        %
        % Args:
        %   name (char): Name of the extra
        %
        % Returns:
        %   logical: Whether the extra value is known and true
            b = isfield(self.extra, name) && self.extra.(name);
        end
        
        function b = isExtraFalse(self, name)
        % Returns whether the value of the given extra is known AND false
        %
        % Args:
        %   name (char): Name of the extra
        %
        % Returns:
        %   logical: Whether the extra value is known and false
            b = isfield(self.extra, name) && ~self.extra.(name);
        end
        
        function b = isExtraSet(self, name)
        % Returns whether the value of the given extra is known
        %
        % Args:
        %   name (char): Name of the extra
        %
        % Returns:
        %   logical: Whether the extra value is known
            b = isfield(self.extra, name);
        end
        
        function U = U(self)
        % Returns the basis of this subrepresentation in its parent
        %
        % Returns:
        %   dense double matrix: The change of basis matrix
            if self.hasCorrection
                U = full(self.D0 * self.U0);
            else
                U = full(self.U0);
            end
        end
        
        function sub1 = nice(self)
        % Tries to recover a nice basis for this subrepresentation
        %
        % Returns:
        %   replab.SubRep: A subrepresentation in the same invariant subspace, 
        %                  but with a rational change of basis if it can be found
            sub1 = replab.rep.niceRep(self);
            if isempty(sub1)
                sub1 = self; % fallback
            end
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % Returns:
        %   dense double matrix: The projector on this invariant subspace,
        %                        expressed in the parent representation
            P = full(self.U0'*(self.D0'*self.D0)*self.U0);
        end
        
        function newSub = collapseParent(self)
        % Returns this subrepresentation as expressed in self.parent.parent
            assert(isa(self.parent, 'replab.SubRep'));
            newU0 = self.U0 * self.parent.U;
            newSub = self.parent.parent.subRepUnitary(newU0, self.extra);
        end
        
        %% Str methods
        
        function s = headerStr(self)
            s = 'Subrepresentation';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.hasCorrection
                for i = 1:size(self.U0, 1)
                    v = replab.str.Normalized(self.U0(i, :), self.D0(i,i));
                    names{1, end+1} = sprintf('U(%d,:)', i);
                    values{1, end+1} = v;
                end
            else
                names{1, end+1} = 'U''';
                values{1, end+1} = self.U';
            end
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Rep(self), {'U0' 'D0'} ...
                );
        end
        
        %% Rep methods
        
        function rho = image(self, g)
            rho = full(self.D0*(self.U0*self.parent.image(g)* self.U0')*self.D0');
        end
        
        function rho = inverseImage(self, g)
            rho = full(self.D0*(self.U0*self.parent.inverseImage(g)* self.U0')*self.D0');
        end
        
    end

end
