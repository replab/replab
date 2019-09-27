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
            error('Abstract');
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
        
    end
    
end
