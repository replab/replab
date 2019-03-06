classdef SdpRep < replab.Str
% A SDP matrix satisfying some symmetry constraints
    
    properties (SetAccess = protected)
        U; % Unitary operator block-diagonalizing the matrix
        vars; % The sdp variables used in the matrix
        nComponents; % The block structure : dimension1 x multiplicity
        multiplicities;
        dimensions1;
        dim; % matrix dimension
    end

    methods
        
        function self = SdpRep(U, dimensions1, multiplicities)
            try
                yalmip('version');
            catch
                error('Yalmip not found in the path');
            end
            self.U = U;
            self.vars = sdpvar(1, sum(multiplicities));
            self.nComponents = length(dimensions1);
            assert(self.nComponents == length(multiplicities));
            self.multiplicities = multiplicities;
            self.dimensions1 = dimensions1;
            self.dim = sum(self.multiplicities.*self.dimensions1);
            assert(self.dim == size(U,1));
        end

    end
    
    methods (Static) % Factory methods
        
        function R = fromGenerators(generators)
            assert(iscell(generators), 'Please specify generators in cell array.');
            
            n = size(generators{1},2);
            group = replab.Permutations(n).subgroup(generators);
            
            irrDecomp = group.naturalRepresentation.irreducible;
            U = irrDecomp.U;
            
            dimensions1 = zeros(1, irrDecomp.nComponents);
            multiplicities = zeros(1, irrDecomp.nComponents);
            for i = 1:irrDecomp.nComponents
                dimensions1(i) = irrDecomp.component(i).dimension1;
                multiplicities(i) = irrDecomp.component(i).multiplicity;
            end
            
            R = replab.SdpRep(U, dimensions1, multiplicities);
        end
        
    end
    
    methods
        
        function M = block(self, i)
        % Returns the desired block
            if (i < 1) || (i > self.nComponents)
                error('Block number out of bound');
            end
            
            
        end
        
        function M = fullMatrix(self)
        % Constructs the full SDP matrix in the natural basis
            d = sdpvar(1);
            d(self.dim) = 0;
            co1 = 0;
            co2 = 0;
            for i = 1:self.nComponents
                d(co2 + [1:self.multiplicities(i)*self.dimensions1(i)]) = ...
                    kron(self.vars(co1+[1:self.multiplicities(i)]), ones(1, self.dimensions1(i)));
                co1 = co1 + self.multiplicities(i);
                co2 = co2 + self.multiplicities(i)*self.dimensions1(i);
            end
            M = self.U*diag(d)*self.U';
        end
        
        function varargout = subsref(self, varargin)
        % Extract one part of the SDP matrix
            switch varargin{1}(1).type
                case '()'
                    [varargout{1:nargout}] = subsref(self.fullMatrix, varargin{1});
                case '.'
                    % This was actually a function call
                    [varargout{1:nargout}] = builtin('subsref',self,varargin{1});
                otherwise
                    error('Not a valid indexing expression');
            end
        end
        
        function see(self)
        % Displays internal info about the matrix composition. For now, we
        % delegate to yalmip
            see(self.fullMatrix);
        end
        
        function val = value(self)
        % Returns the current numerical value of the object
            val = value(self.fullMatrix);
        end
        
%         function y = ge(X,Y)
%         % greater or equal constraint
%             try
%                 y = constraint(self, '>=', Y);
%             catch
%                 error(lasterr)
%             end
%         end
        
        function s = str(self)
        % Nice string representation
            s = ['SDP matrix of size ', num2str(self.dim), 'x', num2str(self.dim), ' with ', num2str(length(self.vars)), ' variables.'];
            s = [s, char(10)];
            s = [s, 'Block structure: '];
            for i = 1:self.nComponents
                s = [s, num2str(self.multiplicities(i)), 'x', num2str(self.dimensions1(i)), ' + '];
            end
            s = s(1:end-3);
        end
        
    end
    
end
