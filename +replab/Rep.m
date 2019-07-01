classdef Rep < replab.Str
% Describes a orthogonal (real) or unitary (complex) finite dimensional represenation
% of a compact group
% 
% Only "image" needs to be implemented in principle; for optimization purposes,
% actions can also be specialized.
    
    properties (SetAccess = protected)
        group;     % Group being represented
        field;     % 'R', 'C' or 'H'
        dimension; % Representation dimension
    end
    
    properties (Access = protected)
        commutant_ = [];
        % irreducible_ = [];
    end
        
    methods % ABSTRACT
        
        function rho = image(self, g) % ABSTRACT
        % Returns the image of the group element
            rho = self.imageFun(g);
        end
        
    end
    
    methods
        
        function c = commutant(self)
        % Returns the commutant of this representation
        % This is the algebra of matrices that commute with
        % this representation images, i.e. for any g in G, we have
        % rho(g) * X = X * rho(g)
            if isequal(self.commutant_, [])
                self.commutant_ = replab.Commutant(self);
            end
            c = self.commutant_;
        end

        % Str

        function s = headerStr(self)
            f = replab.str.field(self.field, 'Orthogonal real', 'Unitary complex', 'Unitary quaternion');
            s = sprintf('%s representation of dimension %d', f, self.dimension);
        end

        % SAMPLING
        
        function rho = sample(self)
            rho = self.image(self.group.sample);
        end
        
        % ACTIONS
        
        function v1 = action(self, g, v)
        % Computes the action of the group element g on the 
        % column vector v through the current representation
        %
        % Note: v can be as well a matrix in which case
        % its columns are treated as independent vectors so that
        % this method computes the matrix(rep)-matrix(element) product
            gI = self.image(g);
            v1 = gI * v;
        end
        
        function X1 = adjointAction(self, g, X)
        % Computes the adjoint action of the group through the current representation
        % For a group element g and a matrix X, it returns rho(g) * X * rho(g)'
            gI = self.image(g);
            X1 = gI * X * gI';
        end
        
        % REAL OR COMPLEX
        
        function complexRep = complexification(self)
            assert(isequal(self.field, 'R'), 'Representation should be real to start with');
            complexRep = replab.RepFun(self.group, 'C', self.dimension, @(g) self.image(g));
        end
        
        function quaternionRep = quaternionification(self)
            assert(isequal(self.field, 'R'), 'Representation should be real to start with');
            quaternionRep = replab.RepFun(self.group, 'H', self.dimension, @(g) self.image(g));
        end
        
        function sub = subRep(self, U)
        % Returns a subrepresentation of this representation,
        % given by the column basis vectors in U, such that
        % sub.image(g) = U' * self.image(g) * U
        %
        % U has dimension self.dimension x dim(subRepresentation)
        %
        % Note that U needs to be orthogonal, but the basis vectors
        % do not need to be normalized
            sub = replab.SubRep(self, U);
        end
        
    end
    
end
