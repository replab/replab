classdef CommutantVar < replab.Str
% CommutantVar is a sdpvar class for matrices satisfying symmetry
% constraints.
%
% See also replab.CommutantVar.fromPermutations,
%          replab.CommutantVar.fromSdpMatrix,
%          replab.CommutantVar.fromSymSdpMatrix

% Warning: this object inherits from a handle object, therefore it is also
% a handle object. To copy this object and obtain two identical but
% independent objects, use the 'copy' method

% Current limitations:
% - The sdp matrix produced is currently always square, real, and symmetric
% - Only supports complex and quaternionic representations with dimension 2
%   and 4 respectively.

% For correct class precedence in matlab, rather use this declaration:
% classdef (InferiorClasses = {?sdpvar,?gem,?sgem}) CommutantVar < replab.Str


    properties (SetAccess = protected)
        % The block structure : dimension1 x multiplicity
        U; % Unitary operator block-diagonalizing the matrix
        nComponents; % Number of block components
        dimensions1; % Dimensions of all irreducible representations
        multiplicities; % Multipliticies of all irreducible representations
        types; % Type of all irreducible representations
        dim; % matrix dimension
        blocks; % The sdp blocks corresponding to each irreducible representation
        linearConstraints; % linear constraints imposed on the matrix
    end

    methods


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                          inherited classes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function s = headerStr(self)
        % s = headerStr(self)
        %
        % Header string representing the object
        %
        % Arg:
        %     self: CommutantVar
        % 
        % Returns:
        %     s: string 
        %
        % See also:
        %     replab.CommutantVar.str
        
            s = sprintf('Commutant variable %dx%d (%d blocks, %d scalar variables)', self.dim, self.dim, length(self.blocks), self.nbVars());
        end

        function names = hiddenFields(self)
        % names = hiddenFields(self)
        %
        % Overload of replab.Str.hiddenFields
        %
        % See also:
        %     replab.Str.hiddenFields
        
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'blocks';
            names{1, end+1} = 'nComponents';
            names{1, end+1} = 'linearConstraints';
        end

        function [names, values] = additionalFields(self)
        % [names, values] = additionalFields(self)
        %
        % Overload of replab.Str.additionalFields
        %
        % See also:
        %     replab.Str.additionalFields

            [names, values] = additionalFields@replab.Str(self);

            names{1, end+1} = 'blocks';
            dims = zeros(1,length(self.blocks));
            for i = 1:length(self.blocks)
                dims(i) = size(self.blocks{i},1);
            end
            values{1, end+1} = dims;

            if ~isempty(self.linearConstraints)
                names{1, end+1} = 'linearConstraints';
                values{1, end+1} = self.linearConstraints;
            end
        end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                               constructors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function self = CommutantVar(generators, sdpMatrix, sdpMatrixIsSym)
        % self = CommutantVar(generators, sdpMatrix, sdpMatrixIsSym)
        %
        % Class constructor, not to be used directly.
        %
        % Args:
        %     generators: a list of generators under which the matrix is to
        %         remain unchanged
        %     sdpMatrix: the SDP matrix on which to impose permutation
        %         invariance (should be empty if none)
        %     sdpMatrixIsSym: whether the structure of the provided
        %         sdpMatrix is already invariant under the generators
        %         (should be empty if no sdpMatrix is provided)
        %
        % Results:
        %     self: replab.CommutantVar object
        %
        % See also:
        %     replab.CommutantVar.fromPermutations,
        %     replab.CommutantVar.fromSdpMatrix,
        %     replab.CommutantVar.fromSymSdpMatrix
        
            try
                yalmip('version');
            catch
                error('Yalmip not found in the path');
            end
            
            % We keep track if some big rounding is done
            maxOuterEpsilonFound = 0;
            maxEpsilonFound = 0;
            epsilonWarning = 1e-10;

            % Input checking
            assert(nargin <= 3, 'Not enough arguments.');
            
            % very special case needed to quickly construct a copy of an
            % object
            if isa(generators, 'replab.CommutantVar') && isempty(sdpMatrix) && isempty(sdpMatrixIsSym)
                rhs = generators;
                % Override properties
%                 fns = properties(rhs);
%                 for i=1:length(fns)
%                     R.(fns{i}) = rhs.(fns{i});
%                 end
                % The above is not supported by octave, so we copy all elements
                % by hand...
                self.U = rhs.U;
                self.nComponents = rhs.nComponents;
                self.dimensions1 = rhs.dimensions1;
                self.multiplicities = rhs.multiplicities;
                self.types = rhs.types;
                self.dim = rhs.dim;
                self.blocks = rhs.blocks;
                self.linearConstraints = rhs.linearConstraints;
                return;
            end

            assert(iscell(generators), 'Please specify generators in cell array.');
            n = size(generators{1}, 2);

            if isempty(sdpMatrix)
                assert(isempty(sdpMatrixIsSym), ['No sdpMatrix provided but sdpMatrixIsSym set to ', num2str(sdpMatrixIsSym)]);
            else
                assert(isequal(size(sdpMatrix), [n n]), 'Wrong matrix or group dimension.');
                assert(isequal(sdpMatrixIsSym, false) || isequal(sdpMatrixIsSym, true), 'sdpMatrixIsSym should be a boolean.');
            end
            
            % Representation decomposition
            group = replab.signed.Permutations(n).subgroup(generators);
            irrDecomp = group.naturalRep.decomposition;
            U = zeros(n, 0);
            dimensions1 = zeros(1,irrDecomp.nComponents);
            multiplicities = zeros(1,irrDecomp.nComponents);
            types = '';
            for i = 1:irrDecomp.nComponents
                component = irrDecomp.component(i);
                dimensions1(i) = component.copyDimension;
                multiplicities(i) = component.multiplicity;
                types(i) = component.copy(1).realDivisionAlgebra.shortName;
                for j = 1:component.multiplicity
                    copy = component.copy(j);
                    % correction, as the new RepLAB convention
                    % is to store basis vectors as row vectors
                    U = [U copy.U'];
                end
            end
            
            % We set most class attributes
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

            % Constructing the SDP blocks now
            if isempty(sdpMatrix) || ~sdpMatrixIsSym
                % We construct the SDP blocks from scratch
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
                
                % We keep in memory the constraint imposed by the SDP matrix
                % TODO: if a SDP matrix was provided, eliminate linear
                %       constraints cleanly by applying the generators
                %       to the SDP matrix. This would be much better.
                if isempty(sdpMatrix)
                    self.linearConstraints = [];
                else
                    assert(isnumeric(sdpMatrix) || isa(sdpMatrix, 'sdpvar'), ['Wrong type for sdpMatrix: ', class(sdpMatrix), '.']);
                    self.linearConstraints = (self.fullMatrix == sdpMatrix);
                end
            else
                % We construct the SDP blocks from the provided SDP matrix
                % off-block-diagonal terms should be zero
                
                assert(isnumeric(sdpMatrix) || isa(sdpMatrix, 'sdpvar'), ['Wrong type for sdpMatrix: ', class(sdpMatrix), '.']);
                
                % We compute each block from the provided SDP matrix
                blockMatrix = self.U'*sdpMatrix*self.U;
                
                co = 0;
                for i = 1:self.nComponents
                    d = self.dimensions1(i);
                    switch self.types(i)
                        case 'R'
                            dimBlock = self.multiplicities(i)*d;
                        case 'C'
                            d = d/2;
                            assert(self.dimensions1(i) == 2, ['Dimension ', num2str(self.dimensions1(i)), ' for a complex representation is currently not supported.']);
                            dimBlock = self.multiplicities(i)*2*d;
                        case 'H'
                            d = d/4;
                            assert(self.dimensions1(i) == 4, ['Dimension ', num2str(self.dimensions1(i)), ' for a quaternionic representation is currently not supported.']);
                            dimBlock = self.multiplicities(i)*4*d;
                    end
                    self.blocks{i} = blockMatrix(co + (1:dimBlock), co + (1:dimBlock));
                    
                    % TODO: check and enforce the fine-grained
                    % block structure for the C and H cases.
                    
                    % Three sanity checks:
                    % 1. block-diagonal form
                    shouldBeZero = blockMatrix(co+(1:dimBlock), co+dimBlock+1:end);
                    indices = [0 getvariables(shouldBeZero)];
                    for ind = indices
                        maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                    end
                    
                    % 2. Fine-grained form of real blocks
                    if isequal(self.types(i), 'R') && (self.dimensions1(i) > 1)
                        m = self.multiplicities(i);
                        tmp = reshape(permute(reshape(self.blocks{i}, [d m d m]), [2 1 4 3]), d*m*[1 1]);
                        for j = 1:d-1
                            for k = j:d
                                if j == k
                                    shouldBeZero = tmp((j-1)*m + (1:m), (k-1)*m + (1:m)) - tmp(1:m,1:m);
                                else
                                    shouldBeZero = tmp((j-1)*m + (1:m), (k-1)*m + (1:m));
                                end
                                indices = [0 getvariables(shouldBeZero)];
                                for ind = indices
                                    maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                                end
                            end
                        end
                        
                        % Enforce the structure
                        self.blocks{i} = tmp(1:m, 1:m);
                    end
                    
                    % 3. hermitianity
                    if ~ishermitian(self.blocks{i})
                        shouldBeZero = self.blocks{i} - self.blocks{i}';
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            maxEpsilonFound = max([maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                        self.blocks{i} = self.blocks{i} + self.blocks{i}';
                    end
                    
                    assert(issymmetric(self.blocks{i}));
                    co = co + dimBlock;
                end
                assert(co == size(sdpMatrix,1))
                
                if maxOuterEpsilonFound > epsilonWarning
                    warning(['The provided SDP matrix does not approximately block-diagonalizes: max delta = ', num2str(maxOuterEpsilonFound), char(10), ...
                        'Off-block-diagonal elements have been replaced by zeros.', char(10), ...
                        'It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                if maxEpsilonFound > epsilonWarning
                    warning(['The provided SDP matrix does not block-diagonalizes into symmetric blocks: max delta = ', num2str(maxEpsilonFound), char(10), ...
                        'Blocks have been symmetrized. It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                
                self.linearConstraints = [];
            end
        end

        function R = copy(rhs)
        % R = copy(rhs)
        %
        % Creates an identical but independent copy of rhs: modifying the 
        % copy does not modify the original object. This is useful
        % internally.
        %
        % Arg:
        %     rhs: replab.CommutantVar object
        %
        % Result:
        %     R: replab.CommutantVar object

            % Create a new simple object
            R = replab.CommutantVar(rhs, [], []);
        end

    end

    methods (Static) % Factory methods

        function R = fromPermutations(generators)
        % R = fromPermutations(generators)
        % 
        % Creates a CommutantVar object: a sdpvar matrix that is invariant
        % under joint permutations of its lines and columns by the provided
        % generators.
        %
        % Args:
        %     generators: list of generators under which the matrix is to
        %         remain unchanged
        %
        % Results:
        %     R: CommutantVar object
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]})
        %
        % See also:
        %     replab.CommutantVar.fromSdpMatrix
        %     replab.CommutantVar.fromSymSdpMatrix
        
            R = replab.CommutantVar(generators, [], []);
        end

        function R = fromSdpMatrix(sdpMatrix, generators)
        % R = fromSdpMatrix(sdpMatrix, generators)
        % 
        % Imposes symmetry constraints onto an existing SDP matrix. This
        % creates a CommutantVar matrix which satisfies both:
        %     - the structure encoded into sdpMatrix
        %     - invariance under joint permutations of its lines and
        %       columns by the provided generators.
        %
        % Args:
        %     sdpMatrix: the SDP matrix on which to impose permutation
        %         invariance
        %     generators: a list of generators under which the matrix is to
        %         remain unchanged
        %
        % Results:
        %     R: a CommutantVar object
        %
        % Example:
        %     sdpMatrix = sdpvar(3);
        %     matrix = replab.CommutantVar.fromSdpMatrix(sdpMatrix, {[3 1 2]})
        %
        % See also:
        %     replab.CommutantVar.fromPermutations
        %     replab.CommutantVar.fromSymSdpMatrix

            R = replab.CommutantVar(generators, sdpMatrix, false);
        end

        function R = fromSymSdpMatrix(sdpMatrix, generators)
        % R = fromSymSdpMatrix(sdpMatrix, generators)
        % 
        % Block-diagonalizes an existing SDP matrix which is already
        % invariant under the group generators.
        %
        % Args:
        %     sdpMatrix: the SDP matrix to block diagonlize
        %     generators: a list of generators under which the matrix
        %         remains unchanged
        %
        % Results:
        %     R: a CommutantVar object
        %
        % Example:
        %     x = sdpvar;
        %     y = sdpvar;
        %     sdpMatrix = [x y y; y x y; y y x];
        %     matrix = replab.CommutantVar.fromSymSdpMatrix(sdpMatrix, {[3 1 2]})
        %
        % See also:
        %     replab.CommutantVar.fromPermutations
        %     replab.CommutantVar.fromSdpMatrix

            R = replab.CommutantVar(generators, sdpMatrix, true);
        end
        
    end


    methods

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                       access to class properties
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function s = str(self)
        % s = str(self)
        %
        % Nice string representation
        %
        % Args:
        %     self: CommutantVar object
        %
        % Results:
        %     s: string
        %
        % See also:
        %     replab.CommutantVar.headerStr
        
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

        function [s1, s2] = size(self, d)
        % [s1, s2] = size(self, [d])
        %
        % Returns the matrix size
        %
        % Args:
        %     self: CommutantVar object
        %     d: specific dimension (optional)
        %
        % Results:
        %     s1: The dimension array, or first dimension if s2 is also
        %         requested
        %     s2: The second dimension (optional)
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]})
        %     size(matrix)
        %
        % See also:
        %     size
        
            if (nargin >= 2) && (d ~= 1) && (d ~= 2)
                error('Wrong dimension in gem::size');
            end

            s1 = self.dim*[1 1];

            if nargin >= 2
                s1 = s1(d);
            elseif nargout == 2
                s2 = s1(2);
                s1 = s1(1);
            end
        end

        function M = block(self, i)
        % M = block(self, i)
        %
        % Returns the desired block
        %
        % Args:
        %     self: CommutantVar object
        %     i: block component number
        % 
        % Results:
        %     M: sdpvar block
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]})
        %     matrix.block(2)
        %
        % See also:
        %     replab.CommutantVar.nComponents
        
            M = self.blocks{i};
        end

        function M = blockMask(self)
        % M = blockMask(self)
        % 
        % Returns a 0-1-filled matrix showing the block structure of the
        % matrix in the irreducible basis.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Results:
        %     M: sparse matrix
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     full(matrix.blockMask)
        
            M = sparse(self.dim, self.dim);
            co = 0;
            for i = 1:self.nComponents
                d = self.dimensions1(i);
                switch self.types(i)
                    case 'R'
                    case 'C'
                        d = d/2;
                    case 'H'
                        d = d/4;
                    otherwise
                        error('Unknown type');
                end
                M(co+[1:d*size(self.blocks{i},1)], co+[1:d*size(self.blocks{i},1)]) = kron(ones(size(self.blocks{i})), eye(d));
                co = co + d*size(self.blocks{i},1);
            end
        end

        function M = fullMatrix(self)
        % M = fullMatrix(self)
        %
        % Constructs the full form of the invariant matrix in the natural
        % basis.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     M: sdpvar matrix
        % 
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     see(matrix.fullMatrix)
        
            M = sdpvar(1);
            M(self.dim^2) = 0;
            M = reshape(M, self.dim*[1 1]);
            co = 0;
            for i = 1:self.nComponents
                d = self.dimensions1(i);
                switch self.types(i)
                    case 'R'
                    case 'C'
                        d = d/2;
                    case 'H'
                        d = d/4;
                    otherwise
                        error('Unknown type');
                end
                M(co + (1:d*size(self.blocks{i},1)), co + (1:d*size(self.blocks{i},1))) = kron(self.blocks{i}, eye(d));
                co = co + d*size(self.blocks{i},1);
            end
            M = self.U*M*self.U';
        end

        function vars = getVariables(self)
        % vars = getVariables(self)
        %
        % Returns the Yalmip indices of the SDP variable used by the
        % object. If the object includes linear constraints, variables from
        % the constraints are also counted.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     vars: vector of indices 
        %
        % Example:
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     matrix1.getVariables
        %     matrix2 = replab.CommutantVar.fromSdpMatrix(sdpvar(3), {[2 3 1]})
        %     matrix2.getVariables
        %     x = sdpvar;
        %     sdpMatrix = [1 x x; x 1 x; x x 1];
        %     matrix3 = replab.CommutantVar.fromSymSdpMatrix(sdpMatrix, {[2 3 1]})
        %     matrix3.getVariables
        %
        % See also:
        %     replab.CommutantVar.nbVars
        %     replab.CommutantVar.getBaseMatrix
        %     sdpvar.getvariables

            vars = getvariables(self.blocks{1});
            for i = 2:self.nComponents
                vars = [vars, getvariables(self.blocks{i})];
            end
            vars = unique(vars);
            
            % If we want to include variables from the linear constraints as well:
            if ~isempty(self.linearConstraints)
                vars = unique([vars, getvariables(self.linearConstraints)]);
            end
        end

        function n = nbVars(self)
        % n = nbVars(self)
        %
        % Returns the number of SDP variable used be the object.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     n: number of SDP variables
        %
        % Example:
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     matrix1.nbVars
        %     matrix2 = replab.CommutantVar.fromSdpMatrix(sdpvar(3), {[2 3 1]})
        %     matrix2.nbVars
        %     x = sdpvar;
        %     sdpMatrix = [1 x x; x 1 x; x x 1];
        %     matrix3 = replab.CommutantVar.fromSymSdpMatrix(sdpMatrix, {[2 3 1]})
        %     matrix3.nbVars
        %
        % See also:
        %     replab.CommutantVar.getVariables
        %     replab.CommutantVar.getBaseMatrix
        
            n = length(self.getVariables);
        end

        function basis = getBaseMatrix(self, index)
        % basis = getBaseMatrix(self, index)
        %
        % Returns the coefficients contributing the the SDP variable with
        % given index.
        %
        % Args:
        %     self: CommutantVar object
        %     index: index of the sdp variable (or 0 for the constant term)
        %
        % Returns:
        %     basis: the matrix of coefficients
        %
        % Example:
        %
        % See also:
        %     replab.CommutantVar.getVariables
        %     replab.CommutantVar.see
        %     sdpvar.getBaseMatrix

            % For now, we simply return the coefficients if they are in the
            % fullMatrix, TODO : also take into account the linear
            % constraint...
            basis = getbasematrix(self.fullMatrix, index);
        end

        function see(self)
        % see(self)
        %
        % Displays internal info about the structure of the matrix in full
        % form.
        %
        % Args:
        %     self: CommutantVar object
        %
        % See also:
        %     sdpvar.see
        
            see(self.fullMatrix);
        end

        function val = value(self)
        % val = value(self)
        %
        % Returns the current numerical value of the object.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     val: numerical value taken by the object
        %
        % See also:
        %     sdpvar.value
        
            val = value(self.fullMatrix);
        end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                            utility methods
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function varargout = subsref(self, varargin)
        % varargout = subsref(self, varargout)
        % 
        % Overload of subsref. This function is called when using the
        % syntax '()' to extract one part of a matrix.
        %
        % Args:
        %     self: CommutantVar object
        %     varargin: usual
        %
        % Results:
        %     varargout: depends on the usage
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     matrix(1)
        %
        % See also:
        %     subsref
        
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

        function okLevel = compatibleWith(X, Y)
        % okLevel = compatibleWith(X, Y)
        %
        % Checks whether the block structure of Y is compatible with that
        % of X.
        % 
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     okLevel:
        %         0 if Y is not compatible with X
        %         1 if Y has the block structure of X
        %         2 if Y has the block structure of X and identical blocks
        %             in X correspond to identical blocks in Y
        %
        % Example:
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     matrix2 = replab.CommutantVar.fromPermutations({[2 1 3]})
        %     % both matrices have comparable block structures ...
        %     full(blockMask(matrix1)) 
        %     full(blockMask(matrix2))
        %     % ... but not in the same basis
        %     matrix1.compatibleWith(matrix2)
        %     % The following matrix has the correct structure in the right basis:
        %     M = matrix1.U*(rand(3).*matrix1.blockMask)*matrix1.U'
        %     matrix1.compatibleWith(M);

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track if some rounding is done
            maxOuterEpsilonFound = 0;
            maxEpsilonFound = 0;
            epsilonWarning = 1e-13;

            % basic tests
            if ~isequal(size(X), size(Y))
                okLevel = 0;
                return;
            end

            % We examine each case independently
            okLevel = 2;
            if isa(X, 'replab.CommutantVar') && (isa(Y, 'replab.CommutantVar') || isa(Y, 'sdpvar'))
                % CommutantVar vs CommutantVar/sdpvar

                % We put Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.fullMatrix*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % We check the structure of rotatedY
                co = 0;
                for i = 1:X.nComponents
                    d = X.dimensions1(i);
                    switch X.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end

                    % First we check the gross structure
                    shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + 1 + d*size(X.blocks{i},2):end);
                    indices = [0 getvariables(shouldBeZero)];
                    for ind = indices
                        if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                            okLevel = 0;
                            break;
                        end
                        maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                    end

                    % Now we check the intro-block structure
                    if okLevel == 2
                        % Check full compatibility
                        ideal = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},2)));
                        ideal = kron(ideal, eye(d));
                        shouldBeZero = ideal - rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2)));

                        % we check if the coefficients are negligible
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                okLevel = 1;
                                maxEpsilonFound = 0;
                                break;
                            end
                            maxEpsilonFound = max([maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                    end
                    if okLevel == 1
                        % Check partial compatibility
                        mask = kron(ones(size(X.blocks{i})), eye(d));
                        shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2))).*(1-mask);

                        % we check if the coefficients are negligible
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                okLevel = 0;
                                return;
                            end
                            maxEpsilonFound = max([maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                    end
                    co = co + d*size(X.blocks{i},1);
                end
            elseif isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar')
                % CommutantVar vs sthg

                % We put Y in the block basis of X
                rotatedY = X.U'*Y*X.U;

                % We check the structure of rotatedY
                co = 0;
                for i = 1:X.nComponents
                    d = X.dimensions1(i);
                    switch X.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end

                    % First we check the gross structure
                    shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + 1 + d*size(X.blocks{i},2):end);
                    if max(max(abs(shouldBeZero))) > epsilon
                        okLevel = 0;
                        return;
                    end
                    maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(shouldBeZero)))]);

                    % Now we check the intro-block structure
                    if okLevel == 2
                        % Check full compatibility
                        ideal = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},2)));
                        ideal = kron(ideal, eye(d));
                        shouldBeZero = ideal - rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2)));
                        if max(max(abs(shouldBeZero))) > epsilon
                            okLevel = 1;
                            maxEpsilonFound = 0;
                        end
                        maxEpsilonFound = max([maxEpsilonFound, max(max(abs(shouldBeZero)))]);
                    end
                    if okLevel == 1
                        % Check partial compatibility
                        mask = kron(ones(size(X.blocks{i})), eye(d));
                        shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2))).*(1-mask);
                        if max(max(abs(shouldBeZero))) > epsilon
                            okLevel = 0;
                            return;
                        end
                        maxEpsilonFound = max([maxEpsilonFound, max(max(abs(shouldBeZero)))]);
                    end
                    co = co + d*size(X.blocks{i},1);
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg vs CommutantVar
                okLevel = Y.compatibleWith(X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end


            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if max([maxOuterEpsilonFound, maxEpsilonFound]) > epsilonWarning
                warning(['Block structure mismatch by ', num2str(maxEpsilonFound), ' ignored.']);
            end
        end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                             algebra methods
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function Z = plus(X,Y)
        % Z = plus(X, Y)
        %
        % addition operator
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     plus
        %     replab.CommutantVar.plus

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-13;

            size1 = size(X);
            size2 = size(Y);

            % Check that dimensions are compatible
            if ~isequal(size1, size2)
                error('Incompatible size for matrix substraction');
            end

        	% We examine each case independently
            if isa(X, 'replab.CommutantVar')
                % CommutantVar + sthg

                % We verify that both variables have fully compatible structures
                compatLevel = X.compatibleWith(Y);
                if compatLevel ~= 2
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                end

                % We express Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.fullMatrix*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % The block structure matches fully, we procede to perform the
                % addition on each block
                Z = copy(X);
                co = 0;
                for i = 1:Z.nComponents
                    d = Z.dimensions1(i);
                    switch Z.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end
                    constantBlock = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},1)));
                    % We make sure that the constant block is hermitian
                    shouldBeZero = constantBlock - constantBlock';
                    if isa(shouldBeZero, 'sdpvar')
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                    else
                        if max(max(abs(shouldBeZero))) > epsilon
                            error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                        end
                        maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                    end
                    if ~ishermitian(constantBlock)
                        % We force exact hermiticity
                        constantBlock = (constantBlock + constantBlock')/2;
                    end
                    Z.blocks{i} = Z.blocks{i} + constantBlock;
                    co = co + d*size(X.blocks{i},1);
                end
                
                % We keep track of the linear constraints
                if isa(Y, 'replab.CommutantVar')
                    Z.linearConstraints = [X.linearConstraints, Y.linearConstraints];
                else
                    Z.linearConstraints = X.linearConstraints;
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg + CommutantVar
                Z = Y+X;
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end

            
            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if maxNonHermiticity > epsilonWarning
                warning(['Non-hermiticity of order ', num2str(maxEpsilonFound), ' was corrected.']);
            end
        end

        function Z = minus(X,Y)
        % Z = minus(X, Y)
        %
        % substraction operator
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     minus
        %     replab.CommutantVar.plus

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-13;

            size1 = size(X);
            size2 = size(Y);

            % Check that dimensions are compatible
            if ~isequal(size1, size2)
                error('Incompatible size for matrix substraction');
            end

        	% We examine each case independently
            if isa(X, 'replab.CommutantVar')
                % CommutantVar - sthg

                % We verify that both variables have fully compatible structures
                compatLevel = X.compatibleWith(Y);
                if compatLevel ~= 2
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                end

                % We express Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.fullMatrix*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % The block structure matches fully, we procede to perform the
                % substraction on each block
                Z = copy(X);
                co = 0;
                for i = 1:Z.nComponents
                    d = Z.dimensions1(i);
                    switch Z.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end
                    constantBlock = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},1)));
                    % We make sure that the constant block is hermitian
                    shouldBeZero = constantBlock - constantBlock';
                    if isa(shouldBeZero, 'sdpvar')
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                    else
                        if max(max(abs(shouldBeZero))) > epsilon
                            error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                        end
                        maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                    end
                    if ~ishermitian(constantBlock)
                        % We force exact hermiticity
                        constantBlock = (constantBlock + constantBlock')/2;
                    end
                    Z.blocks{i} = Z.blocks{i} - constantBlock;
                    co = co + d*size(X.blocks{i},1);
                end
                
                % We keep track of the linear constraints
                if isa(Y, 'replab.CommutantVar')
                    Z.linearConstraints = [X.linearConstraints, Y.linearConstraints];
                else
                    Z.linearConstraints = X.linearConstraints;
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg - CommutantVar
                Z = -(Y-X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end


            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if maxNonHermiticity > epsilonWarning
                warning(['Non-hermiticity of order ', num2str(maxEpsilonFound), ' was corrected.']);
            end
        end

        function X = uminus(self)
        % X = uminus(self)
        %
        % unary minus operator
        %
        % Args:
        %     self: CommutantVar
        %
        % Results:
        %     X: CommutantVar
        %
        % See also:
        %     uminus
        %     replab.CommutantVar.minus

            X = self.copy;
            for i = 1:X.nComponents
                X.blocks{i} = -X.blocks{i};
            end
        end

        function Z = times(X, Y)
        % Z = times(X, Y)
        %
        % Element-wise multiplication. Only for multiplication by a scalar,
        % use fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     times
        %     replab.CommutantVar.rdivide
        %     replab.CommutantVar.fullMatrix

            % Check that dimensions are compatible
            if (numel(X) ~= 1) && (numel(Y) ~= 1)
                error('Use fullMatrix for non-scalar multiplications.');
            end

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % CommutantVar .* CommutantVar
                error('Use fullMatrix instead');
            elseif isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar')
                % CommutantVar .* sthg
                Z = X.copy;
                for i = 1:Z.nComponents
                    Z.blocks{i} = Y.*Z.blocks{i};
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg .* CommutantVar
                Z = Y.*X;
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end
        end

        function Z = rdivide(X, Y)
        % Z = rdivide(X, Y)
        %
        % Element-wise right division operator. Only for division by a
        % scalar, use fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     rdivide
        %     replab.CommutantVar.times
        %     replab.CommutantVar.ldivide
        %     replab.CommutantVar.fullMatrix

            % Check that dimensions are compatible
            if numel(Y) ~= 1
                error('Use fullMatrix for non-scalar multiplications.');
            end
            if ~isa(X, 'replab.CommutantVar') || isa(Y, 'replab.CommutantVar')
                error('Numerator should be a CommutantVar and denominator should not.');
            end

            Z = X.copy;
            for i = 1:Z.nComponents
                Z.blocks{i} = Z.blocks{i}./Y;
            end
        end

        function Z = ldivide(X, Y)
        % Z = ldivide(X, Y)
        %
        % Element-wise left division operator. Only for division by a
        % scalar, use fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     ldivide
        %     replab.CommutantVar.times
        %     replab.CommutantVar.rdivide
        %     replab.CommutantVar.fullMatrix

            Z = Y./X;
        end

        function Z = mtimes(X, Y)
        % Z = mtimes(X, Y)
        %
        % Matrix multiplication. Only for multiplication by a scalar, use
        % fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     mtimes
        %     replab.CommutantVar.mrdivide
        %     replab.CommutantVar.fullMatrix

            % Check that dimensions are compatible
            if (numel(X) ~= 1) && (numel(Y) ~= 1)
                error('Incompatible size for multiplication, use fullMatrix for non-scalar multiplications.');
            end

            Z = X.*Y;
        end

        function Z = mrdivide(X, Y)
        % Z = mrdivide(X, Y)
        %
        % Matrix right division operator. Only for division by a scalar,
        % use fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     mrdivide
        %     replab.CommutantVar.mtimes
        %     replab.CommutantVar.mldivide
        %     replab.CommutantVar.fullMatrix

            % Only scalar division is supported
            if numel(Y) ~= 1
                error('Use fullMatrix for non-scalar division.');
            end

            Z = X./Y;
        end

        function Z = mldivide(X, Y)
        % Z = mldivide(X, Y)
        %
        % Matrix left division operator. Only for division by a scalar,
        % use fullMatrix otherwise.
        %
        % Args:
        %     X: CommutantVar, sdpvar object or numeric matrix
        %     Y: CommutantVar, sdpvar object or numeric matrix
        %
        % Results:
        %     Z: CommutantVar
        %
        % See also:
        %     mldivide
        %     replab.CommutantVar.mtimes
        %     replab.CommutantVar.mrdivide
        %     replab.CommutantVar.fullMatrix

            % Only scalar division is supported
            if numel(X) ~= 1
                error('Use fullMatrix for non-scalar division.');
            end

            Z = Y./X;
        end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                             matrix methods
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function X = trace(self)
        % X = trace(self)
        %
        % Trace operator
        %
        % Args:
        %     self: CommutantVar
        %
        % Results:
        %     X: sdpvar
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     trace(matrix)
        % 
        % See also:
        %     trace

            X = 0;
            for i = 1:self.nComponents
                d = self.dimensions1(i);
                switch self.types(i)
                    case 'R'
                    case 'C'
                        d = d/2;
                    case 'H'
                        d = d/4;
                    otherwise
                        error('Unknown type');
                end
                X = X + trace(self.blocks{i})*d;
            end
        end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                          comparison operators
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function F = eq(X,Y)
        % F = eq(X, Y)
        %
        % Equality constraint. The constraint is imposed through a series
        % of equality constraint for each block. If X or Y includes linear
        % constraints, they are also added to the result.
        %
        % Args:
        %     X: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %     Y: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %
        % Results:
        %     F: Yalmip constraint
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     F = [matrix == 0]
        %     matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(3), {[2 3 1]})
        %     F = [matrix == 0]
        %
        % See also:
        %     eq
        %     replab.CommutantVar.ge

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar') && (numel(Y) == 1)
                % CommutantVar == scalar

                if ~isequal(Y, 0)
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                else
                    F = (X.blocks{1} == Y);
                    for i = 2:X.nComponents
                        F = [F, X.blocks{i} == Y];
                    end
                end
                
                % We add the linear constraints
                F = [F, X.linearConstraints];
            elseif isa(X, 'replab.CommutantVar')
                % CommutantVar == sthg

                % We verify that both variables have compatible structures
                compatLevel = X.compatibleWith(Y);
                if compatLevel == 0
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                end

                % We express Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.U*Y.blockMask*Y.U'*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % We impose the constraint (the constraints might be slightly
                % redundant here, to be improved)
                co = 0;
                F = (sdpvar >= 0);
                for i = 1:X.nComponents
                    d = X.dimensions1(i);
                    switch X.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end
                    for j = 1:1+(2-compatLevel)*(d-1)
                        constantBlock = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},1)));

                        F = [F, X.blocks{i} == constantBlock];

                        if compatLevel == 1
                            co = co + size(X.blocks{i},1);
                        else
                            co = co + d*size(X.blocks{i},1);
                        end
                    end
                end
                F = F(2:end);
                
                % We add the linear constraints
                F = [F, X.linearConstraints];
                if isa(Y, 'replab.CommutatnVar')
                    F = [F, Y.linearConstraints];
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg == CommutantVar
                F = eq(Y,X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end
        end

        function F = ge(X,Y)
        % F = ge(X, Y)
        %
        % Greater or equal semi-definite constraint. The constraint on the
        % matrices are imposed through a series of constraint for each
        % block. If X or Y includes linear constraints, they are also added
        % to the result.
        %
        % Args:
        %     X: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %     Y: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %
        % Results:
        %     F: Yalmip constraint
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     F = [matrix >= 0]
        %     matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(3), {[2 3 1]})
        %     F = [matrix >= 0]
        %
        % See also:
        %     replab.CommutantVar.le

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-13;

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar') && (numel(Y) == 1)
                % CommutantVar >= scalar

                F = (X.blocks{1} >= Y);
                for i = 2:X.nComponents
                    F = [F, X.blocks{i} >= Y];
                end
                
                % We add the linear constraints
                F = [F, X.linearConstraints];
            elseif isa(X, 'replab.CommutantVar')
                % CommutantVar >= sthg

                % We verify that both variables have compatible structures
                compatLevel = X.compatibleWith(Y);
                if compatLevel == 0
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                end

                % We express Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.U*Y.blockMask*Y.U'*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % We impose the constraint (the constraints might be slightly
                % redundant here, to be improved)
                co = 0;
                F = (sdpvar >= 0);
                for i = 1:X.nComponents
                    d = X.dimensions1(i);
                    switch X.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end
                    for j = 1:1+(2-compatLevel)*(d-1)
                        constantBlock = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},1)));

                        % We make sure that the constant block is hermitian
                        shouldBeZero = constantBlock - constantBlock';
                        if isa(shouldBeZero, 'sdpvar')
                            indices = [0 getvariables(shouldBeZero)];
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                                end
                                maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                        end
                        if ~ishermitian(constantBlock)
                            % We force exact hermiticity
                            constantBlock = (constantBlock + constantBlock')/2;
                        end

                        F = [F, X.blocks{i} >= constantBlock];

                        if compatLevel == 1
                            co = co + size(X.blocks{i},1);
                        else
                            co = co + d*size(X.blocks{i},1);
                        end
                    end
                end
                F = F(2:end);
                
                % We add the linear constraints
                F = [F, X.linearConstraints];
                if isa(Y, 'replab.CommutatnVar')
                    F = [F, Y.linearConstraints];
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg >= CommutantVar
                F = le(Y,X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end


            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if maxNonHermiticity > epsilonWarning
                warning(['Non-hermiticity of order ', num2str(maxEpsilonFound), ' was corrected.']);
            end
        end

        function F = le(X,Y)
        % F = le(X, Y)
        %
        % Lesser or equal semi-definite constraint. The constraint on the
        % matrices are imposed through a series of constraint for each
        % block. If X or Y includes linear constraints, they are also added
        % to the result.
        %
        % Args:
        %     X: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %     Y: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %
        % Results:
        %     F: Yalmip constraint
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]})
        %     F = [matrix <= 0]
        %     matrix = replab.CommutantVar.fromSdpMatrix(sdpvar(3), {[2 3 1]})
        %     F = [matrix <= 0]
        %
        % See also:
        %     replab.CommutantVar.ge

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-13;

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar') && (numel(Y) == 1)
                % CommutantVar <= scalar

                F = (X.blocks{1} <= Y);
                for i = 2:X.nComponents
                    F = [F, X.blocks{i} <= Y];
                end
                
                % We add the linear constraints
                F = [F, X.linearConstraints];
            elseif isa(X, 'replab.CommutantVar')
                % CommutantVar <= sthg

                % We verify that both variables have compatible structures
                compatLevel = X.compatibleWith(Y);
                if compatLevel == 0
                    error('Block structure of both matrices don''t match. Consider using fullMatrix.');
                end

                % We express Y in the block basis of X
                if isa(Y, 'replab.CommutantVar')
                    rotatedY = X.U'*Y.U*Y.blockMask*Y.U'*X.U;
                else
                    rotatedY = X.U'*Y*X.U;
                end

                % We impose the constraint (the constraints might be slightly
                % redundant here, to be improved)
                co = 0;
                F = (sdpvar >= 0);
                for i = 1:X.nComponents
                    d = X.dimensions1(i);
                    switch X.types(i)
                        case 'R'
                        case 'C'
                            d = d/2;
                        case 'H'
                            d = d/4;
                        otherwise
                            error('Unknown type');
                    end
                    for j = 1:1+(2-compatLevel)*(d-1)
                        constantBlock = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},1)));

                        % We make sure that the constant block is hermitian
                        shouldBeZero = constantBlock - constantBlock';
                        if isa(shouldBeZero, 'sdpvar')
                            indices = [0 getvariables(shouldBeZero)];
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                                end
                                maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                        end
                        if ~ishermitian(constantBlock)
                            % We force exact hermiticity
                            constantBlock = (constantBlock + constantBlock')/2;
                        end

                        F = [F, X.blocks{i} <= constantBlock];

                        if compatLevel == 1
                            co = co + size(X.blocks{i},1);
                        else
                            co = co + d*size(X.blocks{i},1);
                        end
                    end
                end
                F = F(2:end);

                % We add the linear constraints
                F = [F, X.linearConstraints];
                if isa(Y, 'replab.CommutatnVar')
                    F = [F, Y.linearConstraints];
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg <= CommutantVar
                F = ge(Y,X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end


            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if maxNonHermiticity > epsilonWarning
                warning(['Non-hermiticity of order ', num2str(maxEpsilonFound), ' was corrected.']);
            end
        end

        function F = gt(X,Y)
        % F = gt(X, Y)
        %
        % Greater than constraint. Strict constraints cannot be guarantees,
        % please use ge(X,Y) instead.
        %
        % Args:
        %     X: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %     Y: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %
        % Results:
        %     F: Yalmip constraint
        %
        % See also:
        %     replab.CommutantVar.ge
        
            warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
            F = (X >= Y);
        end

        function F = lt(X,Y)
        % F = lt(X, Y)
        %
        % Lesser than constraint. Strict constraints cannot be guarantees,
        % please use le(X,Y) instead.
        %
        % Args:
        %     X: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %     Y: CommutantVar, compatible sdpvar object or numeric matrix,
        %         or numeric scalar
        %
        % Results:
        %     F: Yalmip constraint
        %
        % See also:
        %     replab.CommutantVar.le
        
            warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
            F = (X <= Y);
        end

    end

end
