classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
%
% The subrepresentation is given by the row vectors of the unitary matrix `U`, such that
% ``self.image(g) = U * parent.image(g) * U'``.

    properties (SetAccess = protected)
        parent % replab.Rep: Parent representation
        U % double matrix, can be sparse: Unitary basis of dimension dChild x dParent
        niceBasis % replab.NiceBasis: Nice decomposition of the basis
        irrepInfo % replab.IrrepInfo or []: Irreducible status information; this representation is known to be
                  %                         irreducible when this field is non empty
    end
    
    methods
        
        function self = SubRep(parent, U, niceBasis, irrepInfo)
        % Constructs a subrepresentation of a parent representation
        %
        % Args:
        %   parent (replab.Rep): Parent representation of which we construct a subrepresentation
        %   U (double matrix, can be sparse): Basis matrix of dimension dThisChild x dParent
        %   niceBasis (replab.NiceBasis or [], optional): Nice decomposition of the basis `U`
        %   irrepInfo (rpelab.IrrepInfo or [], optional): Irreducible status information
            if nargin < 4
                irrepInfo = [];
            end
            if nargin < 3
                niceBasis = [];
            end
            d = size(U, 1);
            dParent = size(U, 2);
            assert(parent.dimension == dParent, 'Incorrect basis dimension');
            assert(isequal(parent.isUnitary, true), 'We support only unitary representations');
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.isUnitary = true;
            self.U = U;
            self.niceBasis = niceBasis;
            self.irrepInfo = irrepInfo;
        end
        
        function b = isKnownIrreducible(self)
        % Returns whether this subrepresentation is known to be irreducible
        %
        % Returns:
        %   logical: True if this subrepresentation is known to be irreducible, 
        %            false if it is reducible or status is unknown
            b = ~isempty(self.irrepInfo);
        end
        
        function b = isKnownCanonicalIrreducible(self)
        % Returns whether this subrepresentation is known to be irreducible and in the canonical division algebra basis
        %
        % Returns:
        %   logical: True if this subrepresentation is known to be irreducible and canonical, 
        %            false if it is reducible/not in the canonical basis, or status is unknown
            if isempty(self.irrepInfo)
                b = false;
                return
            end
            if self.overR
                if isempty(self.irrepInfo.divisionAlgebra)
                    b = false;
                    return
                end
                if isequal(self.irrepInfo.divisionAlgebra, 'R')
                    b = true;
                else
                    b = isequal(self.irrepInfo.isDivisionAlgebraCanonical, true);
                end
            else
                b = true;
                return
            end
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % Returns:
        %   dense double matrix: The projector on this invariant subspace,
        %                        expressed in the parent representation
            P = full(self.U'*self.U);
        end
        
        function newSub = collapseParent(self)
        % Collapses the subrepresentation of a subrepresentation
        %
        % Note that only one level of subrepresentation is removed, as the action is not
        % performed recursively.
        %
        % Returns:
        %   replab.SubRep: A subrepresentation with parent equal to `self.parent.parent`, which
        %                  has in fine the same basis as `self`
            assert(isa(self.parent, 'replab.SubRep'));
            newU = self.U * self.parent.U;
            if ~isempty(self.niceBasis) && ~isempty(self.parent.niceBasis)
                newNiceBasis = self.niceBasis * self.parent.niceBasis;
            else
                newNiceBasis = [];
            end
            newSub = replab.SubRep(self.parent.parent, newU, newNiceBasis, self.irrepInfo);
        end
        
        %% Str methods
        
        function s = headerStr(self)
            if self.isKnownIrreducible
                if self.overR
                    s = 'Real irreducible subrepresentation';
                    if ~isempty(self.irrepInfo.divisionAlgebra)
                        switch self.irrepInfo.divisionAlgebra
                          case 'R'
                            s = [s ' of real type'];
                          case 'C'
                            s = [s ' of complex type'];
                          case 'H'
                            s = [s ' of quaternionic type'];
                        end
                        if self.irrepInfo.isDivisionAlgebraCanonical
                            s = [s ' in the canonical basis'];
                        end
                    end
                else
                    s = 'Complex irreducible subrepresentation';
                end
            else
                s = 'Subrepresentation';
            end
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Rep(self), ...
                {'U'} ...
                );
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if ~isempty(self.niceBasis) && self.niceBasis.isCorrectionDiagonal
                factors = self.niceBasis.normalizationFactors;
                for i = 1:self.dimension
                    v = replab.str.Normalized(self.niceBasis.V(i, :), factors{i});
                    names{1, end+1} = sprintf('U(%d,:)', i);
                    values{1, end+1} = v;
                end
            else
                names{1, end+1} = 'U';
                values{1, end+1} = self.U;
            end
        end
        
        %% Rep methods
        
        function rho = image(self, g)
            rho = full(self.U*self.parent.image(g)*self.U');
        end
        
        function rho = inverseImage(self, g)
            rho = full(self.U*self.parent.inverseImage(g)*self.U');
        end
        
    end

end
