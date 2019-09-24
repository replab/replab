classdef Commutant < replab.Domain
% Describes the algebra of n x n matrices that commute with a group representation
%
% Let rep be a representation of a group G. This describes the set
% of matrices X such that rep.image(g) * X = X * rep.image(g)

    properties (SetAccess = protected)
        parent % replab.Domain: Domain of generic real/complex matrices
        field  % {'R', 'C'}: Vector space defined on real (R) or complex (C) field
        n      % integer: Commutant dimension
        group  % replab.CompactGroup: Group being represented
        rep    % replab.Rep: Representation that this algebra commutes with
    end
    
    properties (Access = protected)
        equivariant_ % replab.Equivariant: Underlying equivariant space
    end
    
    methods (Access = protected)
        
        function e = equivariant(self)            
        % Equivariant space isomorphic to the commutant algebra as a vector space
        %
        % To avoid code duplication, we delegate computations involving the commutant to
        % the equivariant space of a representation to itself.
            if isempty(self.equivariant_)
                self.equivariant_ = replab.makeEquivariant(self.rep, self.rep);
            end
            e = self.equivariant_;
        end
        
    end

    methods

        %% Own methods
        
        function self = Commutant(rep)
            self.rep = rep;
            self.n = rep.dimension;
            self.field = rep.field;
            self.group = rep.group;
            self.parent = replab.domain.Matrices(self.field, self.n, self.n);
        end
        
        function X = project(self, X)
        % Projects any n x n matrix in the invariant subspace
            X = self.equivariant.project(X);
        end

        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d %s commutant matrices', ...
                        self.n, self.n, ...
                        replab.str.field(self.field));
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end
        
        function X = sample(self)
        % Samples a generic matrix from this commutant algebra
            X = self.project(self.parent.sample);
        end
        
        %% Additional sampling methods
        
        function X = sampleSelfAdjoint(self)
        % Samples a generic self-adjoint matrix from this commutant algebra
        %
        % Enforces that X = X'
        % i.e. real matrices are symmetric, complex matrices are Hermitian
            X = self.parent.sample;
            X = (X + X')/2;
            X = self.project(X);
            X = (X + X')/2;
        end
        
        function X = sampleSkewAdjoint(self)
        % Samples a generic skew-adjoint (or anti-self-adjoint) matrix
        %
        % Enforces that X = -X'
            X = self.parent.sample;
            X = (X - X')/2;
            X = self.project(X);
            X = (X - X')/2;
        end

    end
    
end
