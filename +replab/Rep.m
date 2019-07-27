classdef Rep < replab.Str
% Describes a orthogonal (real) or unitary (complex) finite dimensional represenation
% of a compact group
% 
% Only "image" needs to be implemented in principle; for optimization purposes,
% actions can also be specialized.
    
    properties (SetAccess = protected)
        group;     % Group being represented
        field;     % 'R' or 'C'
        dimension; % Representation dimension
    end
    
    properties (Access = protected)
        commutant_ = [];
        decomposition_ = [];
    end
        
    methods % Abstract methods
        
        function rho = image(self, g)
        % Returns the image of the group element
            error('Not implemented');
        end
        
    end
    
    methods % Default implementations
        
        function b = overR(self)
        % Returns true if this representation is defined over the real field
            b = isequal(self.field, 'R');
        end
        
        function b = overC(self)
        % Returns true if this representation is defined over the complex field
            b = isequal(self.field, 'C');
        end
        
        function e = equivariant(self, repC)
        % Returns the space of equivariant linear maps from this rep to repC
        %
        % repR.equivariant(repC) is equivalent to the call 
        % replab.rep.equivariant(repR, repC)
        %
        % the returned space describes the matrices X such that
        %
        % repR.image(g) * X = X * repC.image(g)
            e = replab.rep.equivariant(self, repC)
        end

        function c = commutant(self)
        % Returns the commutant of this representation
        %
        % This is the algebra of matrices that commute with
        % this representation images, i.e. for any g in G, we have
        % rho(g) * X = X * rho(g)
            if isequal(self.commutant_, [])
                self.commutant_ = replab.rep.commutant(self);
            end
            c = self.commutant_;
        end
        
        function I = decomposition(self)
        % Returns the irreducible decomposition of this representation
            if isequal(self.decomposition_, [])
                dec = replab.rep.decomposition(self);
                self.decomposition_ = dec.nice;
            end
            I = self.decomposition_;
        end

        % Str

        function s = headerStr(self)
            f = replab.str.field(self.field, 'Orthogonal real', 'Unitary complex');
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
            complexRep = replab.Rep.lambda(self.group, 'C', self.dimension, @(g) self.image(g));
        end
           
        % Construction of derived representations
        
        function cRep = conj(self)
        % Returns the conjugate representation
            if self.overR
                cRep = self;
            else
                cRep = replab.Rep.lambda(self.group, 'C', self.dimension, @(g) conj(self.image(g)));
            end
        end
        
        function rep = directSum(varargin)
            rep = replab.rep.DirectSumRep(varargin);
        end
        
        function rep = tensor(varargin)
            rep = replab.rep.TensorRep(varargin);
        end
        
        function rep = blkdiag(varargin)
            rep = replab.rep.DirectSumRep(varargin);            
        end
        
        function rep = kron(varargin)
            rep = replab.rep.TensorRep(varargin);
        end
        
        % Manipulation of representation space
        
        function sub = subRep(self, U)
        % Returns a subrepresentation of this representation
        %
        % It is described by the row basis vectors in U, such that
        % sub.image(g) = U' * self.image(g) * U
        %
        % A has dimension self.dimension x dim(subRepresentation)
        %
        % Note that U needs to be orthogonal, but the basis vectors
        % do not need to be normalized
            assert(nargin == 2);
            sub = replab.SubRep(self, U);
        end

        function rep1 = leftConjugate(self, U)
        % Returns the (left) conjugation of this representation
        %
        % A should be a unitary matrix
        %
        % It returns a representation rep1 such that
        %
        % rep1.image(g) = U * self.image(g) * U'
        %
        % Note the opposite convention with respect to subRep
            rep1 = replab.ConjugateRep(U, self);
        end

    end

    methods (Static)
        
        function rep = lambda(group, field, dimension, imageFun)
            rep = replab.lambda.Rep(group, field, dimension, imageFun);
        end
        
    end
    
end
