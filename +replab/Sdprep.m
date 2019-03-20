classdef Sdprep < replab.Str
% A SDP matrix satisfying some symmetry constraints
%
% Current limitations:
% - Signed permutations are not supported yet
% - The sdp matrix produced is currently always square, real, and symmetric
% - Only supports complex and quaternionic representations with dimension 2
%   and 4 respectively.
% - limited arythmetic, essentially we can just:
%   - access the elements of the matrix through the fullMatrix method
%   - define an SDP constraint with a constant bound, such as "M >= 0"

    properties (SetAccess = protected)
        U; % Unitary operator block-diagonalizing the matrix
        nComponents; % The block structure : dimension1 x multiplicity
        dimensions1;
        multiplicities;
        types; % The representation type
        dim; % matrix dimension
        blocks; % The sdp blocks corresponding to each irreducible representation
    end

    methods
        
        function self = Sdprep(U, dimensions1, multiplicities, types)
            try
                yalmip('version');
            catch
                error('Yalmip not found in the path');
            end
            
            % We set the class attributes
            self.U = U;
            self.nComponents = length(dimensions1);
            self.dimensions1 = dimensions1;
            self.multiplicities = multiplicities;
            self.types = types;
            self.dim = sum(self.multiplicities.*self.dimensions1);

            % sanity checks
            assert(self.nComponents == length(self.multiplicities), [num2str(self.nComponents), ' components but ', num2str(length(self.multiplicities)), ' multiplicities']);
            assert(self.nComponents == length(self.types), [num2str(self.nComponents), ' components but ', num2str(length(self.types)), ' types']);
            assert(self.dim == size(self.U,1), ['dimension is ', num2str(self.dim), ' but U is of size ', num2str(size(self.U,1)), 'x', num2str(size(self.U,2))]);
            assert(self.dim == size(self.U,2), ['dimension is ', num2str(self.dim), ' but U is of size ', num2str(size(self.U,1)), 'x', num2str(size(self.U,2))]);

            % We construct the SDP matrix by block
            self.blocks = cell(1,self.nComponents);
            for i = 1:self.nComponents
                switch self.types(i)
                    case 'R'
                        self.blocks{i} = sdpvar(self.multiplicities(i));
                    case 'C'
                        assert(self.dimensions1(i) == 2, ['Dimension ', num2str(self.dimensions1(i)), ' for a complex representation is currently not supported.']);
                        self.blocks{i} = kron(sdpvar(self.multiplicities(i)), eye(2));
                        tmp = kron(sdpvar(self.multiplicities(i), self.multiplicities(i), 'skew'), [1 0; 0 -1]);
                        reOrder = reshape(1:2*self.multiplicities(i), 2, self.multiplicities(i));
                        reOrder = reOrder([2 1],:);
                        reOrder = reOrder(:)';
                        self.blocks{i} = self.blocks{i} + tmp(:, reOrder);
                    case 'H'
                        assert(self.dimensions1(i) == 4, ['Dimension ', num2str(self.dimensions1(i)), ' for a quaternionic representation is currently not supported.']);
                        da = replab.DivisionAlgebra.quaternion;
                        self.blocks{i} = sdpvar(1)*ones(4*self.multiplicities(i));
                        for x = 1:self.multiplicities(i)
                            for y = x:self.multiplicities(i)
                                vars = sdpvar(4,1);
                                if (x == y)
                                    % Only keep the symmetric part
                                    self.blocks{i}(4*(x-1)+[1:4], 4*(y-1)+[1:4]) = da.toMatrix([1 0 0 0]'.*vars);
                                else
                                    % keep everything
                                    self.blocks{i}(4*(x-1)+[1:4], 4*(y-1)+[1:4]) = da.toMatrix(vars);
                                    % Set the off-diagonal terms as well
                                    self.blocks{i}(4*(y-1)+[1:4], 4*(x-1)+[1:4]) = da.toMatrix([1 -1 -1 -1]'.*vars);
                                end
                            end
                        end
                    otherwise
                        error('Unknown type');
                end
                
                % sanity check
                assert(issymmetric(self.blocks{i}));
            end

        end

    end
    
    methods (Static) % Factory methods
        
        function R = fromGenerators(generators)
            assert(iscell(generators), 'Please specify generators in cell array.');
            
            n = size(generators{1},2);
            group = replab.SignedPermutations(n).subgroup(generators);
            
            irrDecomp = group.naturalRepresentation.irreducible;
            U = irrDecomp.U;
            
            dimensions1 = zeros(1, irrDecomp.nComponents);
            multiplicities = zeros(1, irrDecomp.nComponents);
            types = '';
            for i = 1:irrDecomp.nComponents
                dimensions1(i) = irrDecomp.component(i).dimension1;
                multiplicities(i) = irrDecomp.component(i).multiplicity;
                types(i) = irrDecomp.component(i).divisionAlgebra.shortName;
            end
            
            R = replab.Sdprep(U, dimensions1, multiplicities, types);
        end
        
    end
    
    methods
        
        function M = block(self, i)
        % Returns the desired block
            if (i < 1) || (i > self.nComponents)
                error('Block number out of bound');
            end
            
            M = self.blocks{i};
        end
        
        function M = fullMatrix(self)
        % Constructs the full SDP matrix in the natural basis
            M = sdpvar(1);
            M(self.dim^2) = 0;
            M = reshape(M, self.dim*[1 1]);
            co = 0;
            for i = 1:self.nComponents
                dim = self.dimensions1(i);
                switch self.types(i)
                    case 'R'
                    case 'C'
                        dim = dim/2;
                    case 'H'
                        dim = dim/4;
                    otherwise
                        error('Unknown type');
                end
                M(co+[1:dim*size(self.blocks{i},1)], co+[1:dim*size(self.blocks{i},1)]) = kron(self.blocks{i}, eye(dim));
                co = co + dim*size(self.blocks{i},1);
            end
            M = self.U*M*self.U';
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
        
        function F = ge(X,Y)
        % greater or equal constraint
            % We only support some simple cases for now
            if isa(X, 'replab.Sdprep')
                self = X;
                other = Y;
                if isa(Y, 'replab.Sdprep')
                    error('Comparison between two replab.Sdprep not yet supported');
                elseif numel(Y) ~= 1
                    error('Only comparison with scalar supported');
                end
            elseif isa(Y, 'replab.Sdprep')
                self = Y;
                other = X;
                if numel(X) ~= 1
                    error('Only comparison with scalar supported');
                end
            else
                error('Neither of the two arguments is of type Sdprep');
            end
            
            F = [self.blocks{1} >= other];
            for i = 2:self.nComponents
                F = [F, self.blocks{i} >= other];
            end
        end
        
        function s = str(self)
        % Nice string representation
            s = ['SDP matrix of size ', num2str(self.dim), 'x', num2str(self.dim), ' with ', num2str(self.nbVars), ' variables.'];
            s = [s, char(10)];
            s = [s, 'Block structure: '];
            for i = 1:self.nComponents
                switch self.types(i)
                    case 'R'
                        s = [s, num2str(self.dimensions1(i)), '*', num2str(self.multiplicities(i)), 'x', num2str(self.multiplicities(i)), ' + '];
                    case 'C'
                        s = [s, num2str(self.dimensions1(i)/2), '*', num2str(2*self.multiplicities(i)), 'x', num2str(2*self.multiplicities(i)), ' + '];
                    case 'H'
                        s = [s, num2str(self.dimensions1(i)/4), '*', num2str(4*self.multiplicities(i)), 'x', num2str(4*self.multiplicities(i)), ' + '];
                    otherwise
                        error('Unknown type');
                end
            end
            s = s(1:end-3);
        end
        
        function n = nbVars(self)
        % Returns the number of SDP variable used be the object
            n = 0;
            for i = 1:self.nComponents
                n = n + length(getvariables(self.blocks{i}));
            end
        end
        
    end
    
    
end
