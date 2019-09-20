classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
%
% The basis of this subrepresentation in the parent representation is formed by the rows of the
% matrix ``U = diag(D0) * U0``, and we have ``image(g) = U * parent.image(g) * U'``,
% where U0 is an orthogonal basis, not necessarily normalized, and D0 is a correction factor vector.

    properties (SetAccess = protected)
        parent % replab.Rep: Parent representation
        U0 % double matrix, can be sparse: Orthogonal basis of dimension dChild x dParent
        D0 % double row vector: correction factors of dimension 1 x dChild
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
            D0 = ones(1, d);
            hasCorrection = false;
            for i = 1:d
                nrm = self.U0(i,:) * self.U0(i,:)';
                if abs(nrm - 1) > replab.Settings.doubleEigTol
                    hasCorrection = true;
                    D0(i) = 1/sqrt(nrm);
                end
            end
            if hasCorrection
                self.D0 = D0;
            else
                self.D0 = [];
            end
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
        
        function b = hasCorrection(self)
        % Returns true if the basis is not normalized and needs correction
            b = ~isequal(self.D0, []);
        end
        
        function U = U(self)
        % Returns the basis of this subrepresentation in its parent
            if self.hasCorrection
                U = diag(self.D0) * full(self.U0);
            else
                U = full(self.U0);
            end
        end
        
        function sub1 = nice(self)
        % Tries to recover a nice basis for this subrepresentation
            sub1 = replab.rep.niceRep(self);
            if isempty(sub1)
                sub1 = self; % fallback
            end
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % The projector is expressed in the parent representation
            P = self.U'*self.U;
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
                    v = replab.str.Normalized(self.U0(i, :), self.D0(i));
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
            if self.hasCorrection
                D = diag(self.D0);
                rho = D * self.U0 * self.parent.image(g) * self.U0' * D';
            else
                rho = self.U0 * self.parent.image(g) * self.U0';
            end
            rho = full(rho);
        end
        
        function rho = inverseImage(self, g)
            if self.hasCorrection
                D = diag(self.D0);
                rho = D * self.U0 * self.parent.inverseImage(g) * self.U0' * D';
            else
                rho = self.U0 * self.parent.inverseImage(g) * self.U0';
            end
            rho = full(rho);
        end
        
    end

end
