classdef Equivariant < replab.Domain
% Describes a vector space of group-equivariant matrices
%
% Let repR and repC be two representations of the same group G.
%
% This describes the set of matrices X such that repR.image(g) * X = X * repC.image(g)
%
% See Proposition 4 of 
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
    
    properties (SetAccess = protected)
        field % {'R', 'C'}: field of the vector space real (R) or complex x(C)
        nR % integer: row size
        nC % integer: column size
        group % replab.CompactGroup: group being represented
        repR % replab.Rep: representation of row space
        repC % replab.Rep: representation of column space
    end
    
    properties (Access = protected)
        parent_ % parent domain, real or complex matrices
    end
    
    methods

        function X = project(self, X)
        % Projects any nR x nC matrix in the equivariant subspace
            error('Abstract');
        end

        function self = Equivariant(repR, repC)
        % Constructor; please do not call this from user code, but
        % rather use `replab.EquivariantDispatch.instance.call(repR, repC)`, 
        % which can eventually select an optimized implementation 
        % depending on the use case.
            self.repR = repR;
            self.nR = repR.dimension;
            self.repC = repC;
            self.nC = repC.dimension;
            assert(isequal(repR.field, repC.field), ...
                   'Both representations must have be defined on the same field');
            self.field = repR.field;
            assert(repR.group == repC.group, ...
                   'Both representations must be defined on the same group');
            self.group = repR.group;
            switch self.field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(self.nR, self.nC);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(self.nR, self.nC);
              otherwise
                error('Unknown field');
            end
        end

        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d %s equivariant matrices', ...
                        self.nR, self.nC, ...
                        replab.str.field(self.field));
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.project(self.parent_.sample);
        end
        
    end
    
end
