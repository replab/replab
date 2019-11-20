classdef CommutantVar < replab.Str
% CommutantVar is a class of sdpvar matrices satisfying symmetry
% constraints. The matrices produced are always invariant under
% transposition.
% 
% A general matrix only subject to invariance under a permutation group can
% be constructed with the method fromPermutations. fromIndexMatrix allows
% one to construct a symmetry-invariant matrix with additional structure.
% Symmetry constraints can also be imposed on existing sdpvar object with
% the fromSdpMatrix constructor. If the provided SDP matrix is already
% known to be invariant under the group, then fromSymSdpMatrix can be used
% to only add the induced block structure onto this matrix.
% 
% See also replab.CommutantVar.fromPermutations,
%          replab.CommutantVar.fromIndexMatrix,
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
        U_; % Unitary operator block-diagonalizing the matrix
        nComponents; % Number of block components
        dimensions1; % Dimensions of all irreducible representations
        multiplicities; % Multipliticies of all irreducible representations
        types; % Type of all irreducible representations
        dim; % matrix dimension
        blocks; % The sdp blocks corresponding to each irreducible representation
        linearConstraints; % linear constraints imposed on the matrix
        fullBlockMatrix_; % The combinations, lazy evaluated
        matrixType; % full, symmetric of hermitian
        field; % real or complex
        sdpMatrix_; % The provided or constructed sdpMatrix
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
            names{1, end+1} = 'fullBlockMatrix_';
            names{1, end+1} = 'U_';
            names{1, end+1} = 'sdpMatrix_';
        end

        function [names, values] = additionalFields(self)
        % [names, values] = additionalFields(self)
        %
        % Overload of replab.Str.additionalFields
        %
        % See also:
        %     replab.Str.additionalFields

            [names, values] = additionalFields@replab.Str(self);

            names{1,end+1} = 'U';
            values{1,end+1} = self.U;

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

        function self = CommutantVar(generators, sdpMatrix, sdpMatrixIsSym, matrixType, field)
        % self = CommutantVar(generators, sdpMatrix, sdpMatrixIsSym, matrixType, field)
        %
        % Class constructor, not to be used directly.
        %
        % Args:
        %     generators: a list of generators under which the matrix is to
        %         remain unchanged
        %     sdpMatrix: the SDP matrix on which to impose permutation
        %         invariance (should be empty if none)
        %     sdpMatrixIsSym: (should be empty if no sdpMatrix is provided)
        %         0: if the sdpMatrix is know expected to be invariant
        %             under the group action
        %         1: if the sdpMatrix is expected to be invariant under the
        %             group action but we wish to check it
        %         2: if the sdpMatrix is known to be invariant under the
        %             group action. In this case, no checks are made
        %     matrixType: one of the following:
        %         'full' : no particular structure
        %         'symmetric' : transpose-invariant matrix
        %         'hermitian' : hermitian matrix
        %     field: matrix elements are either 'real' or 'complex'
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
            maxOuterEpsilonFound = 0; % Out of the block structure elements
            maxEpsilonFoundCH = 0; % If the internal complex/quaternionic structure is not exact
            maxEpsilonFoundField = 0; % If real numbers are not exactly real
            maxEpsilonFoundType = 0; % If symmetric/hermitian matrix structure is not exact
            epsilonWarning = 1e-10;

            % very special case needed to quickly construct a copy of an
            % object
            if isa(generators, 'replab.CommutantVar') && isempty(sdpMatrix) && isempty(sdpMatrixIsSym) && isempty(matrixType) && isempty(field)
                rhs = generators;
                % Override properties
%                 fns = properties(rhs);
%                 for i=1:length(fns)
%                     R.(fns{i}) = rhs.(fns{i});
%                 end
                % The above is not supported by octave, so we copy all elements
                % by hand...
                self.U_ = rhs.U_;
                self.nComponents = rhs.nComponents;
                self.dimensions1 = rhs.dimensions1;
                self.multiplicities = rhs.multiplicities;
                self.types = rhs.types;
                self.dim = rhs.dim;
                self.blocks = rhs.blocks;
                self.linearConstraints = rhs.linearConstraints;
                self.matrixType = rhs.matrixType;
                self.field = rhs.field;
                self.sdpMatrix_ = rhs.sdpMatrix_;
                return;
            end

            % Input checking
            assert(nargin <= 5, 'Not enough arguments.');
            
            assert(iscell(generators), 'Please specify generators in cell array.');
            n = size(generators{1}, 2);

            if isempty(sdpMatrix)
                assert(isempty(sdpMatrixIsSym), ['No sdpMatrix provided but sdpMatrixIsSym set to ', num2str(sdpMatrixIsSym)]);
            else
                assert(isequal(size(sdpMatrix), [n n]), 'Wrong matrix or group dimension.');
                assert(isequal(sdpMatrixIsSym, 0) || isequal(sdpMatrixIsSym, 1) || isequal(sdpMatrixIsSym, 2), 'sdpMatrixIsSym can only take value 0,1,2.');
            end
            
            assert(isequal(matrixType, 'full') || isequal(matrixType, 'symmetric') || isequal(matrixType, 'hermitian'), 'The matrix type must be ''full'', ''symmetric'' or ''hermitian''.');
            assert(isequal(field, 'real') || isequal(field, 'complex'), 'The field must be ''real'' or ''complex''.');
            if isequal(matrixType, 'hermitian') && isequal(field, 'real')
                matrixType = 'symmetric';
            end
            self.matrixType = matrixType;
            self.field = field;
            
            % Representation decomposition
            group = replab.signed.Permutations(n).subgroup(generators);
            irrDecomp = group.definingRep.decomposition;
            U = zeros(n, 0);
            dimensions1 = zeros(1,irrDecomp.nComponents);
            multiplicities = zeros(1,irrDecomp.nComponents);
            types = '';
            for i = 1:irrDecomp.nComponents
                component = irrDecomp.component(i);
                dimensions1(i) = component.irrepDimension;
                multiplicities(i) = component.multiplicity;
                types(i) = component.irrep(1).irrepInfo.divisionAlgebra;
                for j = 1:component.multiplicity
                    copy = component.irrep(j);
                    % correction, as the new RepLAB convention
                    % is to store basis vectors as row vectors
                    U = [U copy.U'];
                end
            end
            
            % We set most class attributes
            self.U_ = U;
            self.nComponents = length(dimensions1);
            self.dimensions1 = dimensions1;
            self.multiplicities = multiplicities;
            self.types = types;
            self.dim = sum(self.multiplicities.*self.dimensions1);

            % sanity checks
            assert(self.nComponents == length(self.multiplicities), [num2str(self.nComponents), ' components but ', num2str(length(self.multiplicities)), ' multiplicities']);
            assert(self.nComponents == length(self.types), [num2str(self.nComponents), ' components but ', num2str(length(self.types)), ' types']);
            assert(self.dim == size(self.U_,1), ['dimension is ', num2str(self.dim), ' but U is of size ', num2str(size(self.U_,1)), 'x', num2str(size(self.U_,2))]);
            assert(self.dim == size(self.U_,2), ['dimension is ', num2str(self.dim), ' but U is of size ', num2str(size(self.U_,1)), 'x', num2str(size(self.U_,2))]);

            % Constructing the SDP blocks now
            if isempty(sdpMatrix) || (sdpMatrixIsSym == 0)
                % We construct the SDP blocks from scratch
                self.blocks = cell(1,self.nComponents);
                for i = 1:self.nComponents
                    switch self.types(i)
                        case 'R'
                            if isequal(matrixType, 'full') && isequal(field, 'real')
                                self.blocks{i} = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'real')
                                self.blocks{i} = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'real');
                            elseif isequal(matrixType, 'full') && isequal(field, 'complex')
                                self.blocks{i} = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'complex')
                                self.blocks{i} = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'complex');
                            elseif isequal(matrixType, 'hermitian') && isequal(field, 'complex')
                                self.blocks{i} = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'complex');
                            end
                            
                        case 'C'
                            basis1 = [1  0   % from replab.domain.ComplexTypeMatrices.toMatrix(1,0);
                                      0  1];
                            basis2 = [0 -1   % from replab.domain.ComplexTypeMatrices.toMatrix(0,1);
                                      1  0];
                            
                            if isequal(matrixType, 'full') && isequal(field, 'real')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'real')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'real');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars2 = imag(vars2); % this part should be antisymmetric
                            elseif isequal(matrixType, 'full') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'complex');
                                vars2R = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'complex');
                                vars2I = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'complex');
                                vars2 = imag(vars2R) + 1i*imag(vars2I); % this part should be fully antisymmetric
                            elseif isequal(matrixType, 'hermitian') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'complex');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'complex');
                                vars2 = 1i*vars2; % Real part should be antisymmetric, and imaginary part symmetric
                            end
                            
                            self.blocks{i} = kron(vars1, basis1) + kron(vars2, basis2);
                            
                        case 'H'
                            basis1 = [1  0  0  0   % from replab.domain.QuaternionTypeMatrices.toMatrix(1,0,0,0);
                                      0  1  0  0
                                      0  0  1  0
                                      0  0  0  1];
                            basis2 = [0 -1  0  0   % from replab.domain.QuaternionTypeMatrices.toMatrix(0,1,0,0);
                                      1  0  0  0
                                      0  0  0 -1
                                      0  0  1  0];
                            basis3 = [0  0 -1  0   % from replab.domain.QuaternionTypeMatrices.toMatrix(0,0,1,0);
                                      0  0  0  1
                                      1  0  0  0
                                      0 -1  0  0];
                            basis4 = [0  0  0 -1   % from replab.domain.QuaternionTypeMatrices.toMatrix(0,0,0,1);
                                      0  0 -1  0
                                      0  1  0  0
                                      1  0  0  0];
                            
                            if isequal(matrixType, 'full') && isequal(field, 'real')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                                vars3 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                                vars4 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'real');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'real')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'real');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars2 = imag(vars2); % this part should be antisymmetric
                                vars3 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars3 = imag(vars3); % this part should be antisymmetric
                                vars4 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars4 = imag(vars4); % this part should be antisymmetric
                            elseif isequal(matrixType, 'full') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                                vars3 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                                vars4 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'full', 'complex');
                            elseif isequal(matrixType, 'symmetric') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'real');
                                vars2R = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars2I = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars2 = imag(vars2R) + 1i*imag(vars2I); % this part should be fully antisymmetric
                                vars3R = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars3I = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars3 = imag(vars3R) + 1i*imag(vars3I); % this part should be fully antisymmetric
                                vars4R = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars4I = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars4 = imag(vars4R) + 1i*imag(vars4I); % this part should be fully antisymmetric
                            elseif isequal(matrixType, 'hermitian') && isequal(field, 'complex')
                                vars1 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'symmetric', 'real');
                                vars2 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars2 = 1i*vars2; % Real part should be antisymmetric, and imaginary part symmetric
                                vars3 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars3 = 1i*vars3; % Real part should be antisymmetric, and imaginary part symmetric
                                vars4 = sdpvar(self.multiplicities(i), self.multiplicities(i), 'hermitian', 'real');
                                vars4 = 1i*vars4; % Real part should be antisymmetric, and imaginary part symmetric
                            end
                            
                            self.blocks{i} = kron(vars1, basis1) + kron(vars2, basis2) + kron(vars3, basis3) + kron(vars4, basis4);

                        otherwise
                            error('Unknown type');
                    end
                end
                
                % We keep in memory the constraint imposed by the SDP matrix
                % TODO: if a SDP matrix was provided, eliminate linear
                %       constraints cleanly by applying the generators
                %       to the SDP matrix. This would be much better.
                if isempty(sdpMatrix)
                    self.sdpMatrix_ = [];
                    self.linearConstraints = [];
                else
                    assert(isnumeric(sdpMatrix) || isa(sdpMatrix, 'sdpvar'), ['Wrong type for sdpMatrix: ', class(sdpMatrix), '.']);
                    self.sdpMatrix_ = sdpMatrix;
                    self.linearConstraints = (self.U_*self.fullBlockMatrix*self.U_' == sdpMatrix);
                end
            else
                % We construct the SDP blocks from the provided SDP matrix
                % off-block-diagonal terms should be zero but will be
                % checked
                
                assert(isnumeric(sdpMatrix) || isa(sdpMatrix, 'sdpvar'), ['Wrong type for sdpMatrix: ', class(sdpMatrix), '.']);
                
                % We compute each block from the provided SDP matrix
                co = 0;
                for i = 1:self.nComponents
                    d = self.dimensions1(i);
                    m = self.multiplicities(i);
                    dimBlock = m*d;
                    self.blocks{i} = self.U_(:, co + (1:dimBlock))' * sdpMatrix * self.U_(:, co + (1:dimBlock));
                    
                    % We also compute the structure when the field
                    % representation dimension is associated with the
                    % multiplicity
                    switch self.types(i)
                        case 'R'
                            mp = m;
                            dp = d;
                        case 'C'
                            mp = m*2;
                            dp = d/2;
                        case 'H'
                            mp = m*4;
                            dp = d/4;
                    end

                    % Five sanity checks (plus reducing the block size if possible):
                    % 1. block-diagonal form
                    if (sdpMatrixIsSym == 1)
                        shouldBeZero = (self.U_(:,co+(1:dimBlock))' * sdpMatrix) * self.U_(:,co+dimBlock+1:end);
                        indices = [0 getvariables(shouldBeZero)];
                        for ind = indices
                            maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                    end
                    
                    % 2. Fine-grained form of blocks in presence of large
                    % dimensions
                    if (dp > 1)
                        tmp = reshape(permute(reshape(self.blocks{i}, [dp mp dp mp]), [2 1 4 3]), dp*mp*[1 1]);
                        
                        if (sdpMatrixIsSym == 1)
                            for j = 1:dp-1
                                for k = j:dp
                                    if j == k
                                        shouldBeZero = tmp((j-1)*mp + (1:mp), (k-1)*mp + (1:mp)) - tmp(1:mp,1:mp);
                                    else
                                        shouldBeZero = tmp((j-1)*mp + (1:mp), (k-1)*mp + (1:mp));
                                    end
                                    indices = [0 getvariables(shouldBeZero)];
                                    for ind = indices
                                        maxOuterEpsilonFound = max([maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                                    end
                                end
                            end
                        end
                        
                        % Enforce the structure
                        self.blocks{i} = tmp(1:mp, 1:mp);
                    end
                    
                    % 3. Check the block structure for complex and
                    % quaternionic representations
                    if (sdpMatrixIsSym == 1) && ~isequal(self.types(i), 'R')
                        if isequal(self.types(i), 'C')
                            checkbases = {[1 0; 0 -1], [0 1; 1 0]};
                        elseif isequal(self.types(i), 'H')
                            checkbases = {[1 0 0 0; 0 -1 0 0; 0 0 0 0], [1 0 0 0; 0 0 0 0; 0 0 -1 0; 0 0 0 0], [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 -1], ...
                                [0 1 0 0; 1 0 0 0; 0 0 0 0; 0 0 0 0], [0 1 0 0; 0 0 0 0; 0 0 0 -1; 0 0 0 0], [0 1 0 0; 0 0 0 0; 0 0 0 0; 0 0 1 0], ...
                                [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0], [0 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0], [0 0 1 0; 0 0 0 0; 0 0 0 0; 0 -1 0 0], ...
                                [0 0 0 1; 0 0 -1 0; 0 0 0 0; 0 0 0 0], [0 0 0 1; 0 0 0 0; 0 1 0 0; 0 0 0 0], [0 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 0]};
                        end
                        for j = 1:length(checkbases)
                            % We just need to check the first element
                            shouldBeZero = sum(sum( checkbases{j}.*self.blocks{i}(1:size(checkbases{j},1), 1:size(checkbases{j},2)) ));
                            indices = [0 getvariables(shouldBeZero)];
                            for ind = indices
                                maxEpsilonFoundCH = max([maxEpsilonFoundCH, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                            end
                        end
                    end
                
                    % 4. Check for complexity
                    if (sdpMatrixIsSym == 1) && isequal(field, 'real')
                        shouldBeZero = imag(self.blocks{i});
                        indices = [0 getvariables(shouldBeZero)];
                        tmp = 0;
                        for ind = indices
                            tmp = max([tmp, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                        maxEpsilonFoundType = max([maxEpsilonFoundType, tmp]);
                        if tmp > 0
                            self.blocks{i} = real(self.blocks{i});
                        end
                    end
                    
                    % 5. Check symmetry or hermiticity
                    if (sdpMatrixIsSym == 1) && ~isequal(matrixType, 'full')
                        if isequal(matrixType, 'symmetric')
                            shouldBeZero = self.blocks{i} - self.blocks{i}.';
                        elseif isequal(matrixType, 'hermitian')
                            shouldBeZero = self.blocks{i} - self.blocks{i}';
                        end
                        indices = [0 getvariables(shouldBeZero)];
                        tmp = 0;
                        for ind = indices
                            tmp = max([tmp, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                        end
                        maxEpsilonFoundType = max([maxEpsilonFoundType, tmp]);
                        if tmp > 0
                            if isequal(matrixType, 'symmetric')
                                self.blocks{i} = (self.blocks{i} + self.blocks{i}.')/2;
                            elseif isequal(matrixType, 'hermitian')
                                self.blocks{i} = (self.blocks{i} + self.blocks{i}')/2;
                            end
                        end
                        if isequal(matrixType, 'symmetric')
                            assert(issymmetric(self.blocks{i}) == 1);
                        elseif isequal(matrixType, 'hermitian')
                            assert(ishermitian(self.blocks{i}) == 1);
                        end
                    end
                    
                    co = co + dimBlock;
                end
                assert(co == size(sdpMatrix,1))
                
                if maxOuterEpsilonFound > epsilonWarning
                    warning(['The provided SDP matrix does not approximately block-diagonalizes: max delta = ', num2str(maxOuterEpsilonFound), char(10), ...
                        'Off-block-diagonal elements have been replaced by zeros.', char(10), ...
                        'It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                if maxEpsilonFoundCH > epsilonWarning
                    warning(['The provided SDP matrix does not block-diagonalizes into blocks of the expected real/complex/quaternionic structure: max delta = ', num2str(maxEpsilonFoundCH), char(10), ...
                        'The structure has not been enforced.', char(10), ...
                        'It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                if maxEpsilonFoundField > epsilonWarning
                    warning(['The provided SDP matrix does not block-diagonalizes into real blocks: max delta = ', num2str(maxEpsilonFoundField), char(10), ...
                        'Imaginary part has been set to zero.', char(10), ...
                        'It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                if maxEpsilonFoundType > epsilonWarning
                    warning(['The provided SDP matrix does not block-diagonalizes into symmetric/hermitian form: max delta = ', num2str(maxEpsilonFoundType), char(10), ...
                        'Blocks have been symmetrized.', char(10), ...
                        'It might be better to use replab.CommutantVar.fromSdpMatrix.']);
                end
                
                self.sdpMatrix_ = sdpMatrix;
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
            R = replab.CommutantVar(rhs, [], [], [], []);
        end

    end

    methods (Static) % Factory methods

        function R = fromPermutations(generators, matrixType, field)
        % R = fromPermutations(generators, matrixType, field)
        % 
        % Creates a CommutantVar object: a sdpvar matrix that is invariant
        % under joint permutations of its lines and columns by the provided
        % generators.
        %
        % Args:
        %     generators: list of generators under which the matrix is to
        %         remain unchanged
        %     matrixType: one of the following:
        %         'full' : no particular structure
        %         'symmetric' : transpose-invariant matrix
        %         'hermitian' : hermitian matrix
        %     field: matrix elements are either 'real' or 'complex'
        %
        % Results:
        %     R: CommutantVar object
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]}, 'symmetric', 'real')
        %
        % See also:
        %     replab.CommutantVar.fromIndexMatrix
        %     replab.CommutantVar.fromSdpMatrix
        %     replab.CommutantVar.fromSymSdpMatrix
        
            R = replab.CommutantVar(generators, [], [], matrixType, field);
        end
        
        function subsets = burn(pairs)
        % subsets = burn(pairs)
        %
        % Performs the burning algorithm on the network described by the
        % edges given in pairs.
        %
        % Args:
        %     pairs: nx2 array of vertices linked by an edge
        %
        % Returns:
        %     subsets: cell array with connex components
        
            uniquesLeft = unique(pairs);
            subsets = {};
            co1 = 0;
            while ~isempty(uniquesLeft)
                co1 = co1 + 1;
                set = uniquesLeft(1);
                co2 = 0;
                while co2 < length(set)
                    co2 = co2 + 1;
                    element = set(co2);
                    sel1 = find(pairs(:,1) == element);
                    sel2 = find(pairs(:,2) == element);
                    newElements = unique([pairs(sel1,2); pairs(sel2,1)])';
                    newElements = setdiff(newElements, set);
                    set = [set, newElements];
                end
                subsets{co1} = sort(set);
                uniquesLeft = setdiff(uniquesLeft, set);
            end
            %subsets{:}
        end
        
        function [ok, coloring] = graphIsBipartite(pairs)
        % [ok, coloring] = graphIsBipartite(pairs)
        %
        % Checks whether the graph described by a list of pairs of
        % connected vertices is bipartite, i.e. whether the graph can be
        % colored with 2 colors.
        %
        % Args:
        %     pairs: nx2 array listing the (undirected) edges connecting
        %         pairs of vertices. Should only contain strictly positive
        %         integers.
        %
        % Returns:
        %     ok: 1 if the graph is bipartite
        %         0 if the graph requires more than 2 colors to be
        %           colored
        %     coloring: coloring of the vertices in case of success
        
            % We list all vertices
            vertices = unique(pairs(:))';
            
            % We assign a neutral color to every vertex
            colors = zeros(1,max(vertices));
            
            coloredVertices = zeros(size(vertices));
            coColoredVertices = 0;
            coV = 0;
            ok = 0;
            while coColoredVertices < length(vertices)
                % We assign color 1 to one new un-reached vertex
                newVertex = vertices(find(setdiff(vertices, coloredVertices), 1, 'first'));
                colors(newVertex) = 1;
                coColoredVertices = coColoredVertices + 1;
                coloredVertices(coColoredVertices) = newVertex;

                % We assign alternating colors to neighbours
                while coV < length(coloredVertices)
                    coV = coV + 1;

                    % We pick the next 
                    currentVertex = coloredVertices(coV); % current vertex
                    currentColor = colors(currentVertex); % color of current vertex

                    % We color all neighbours of this vertex
                    for i = 1:size(pairs,1)
                        newVertex = 0; % new vertex
                        newColor = mod(currentColor,2)+1; % color of neighbouring vertices

                        if pairs(i,1) == currentVertex
                            newVertex = pairs(i,2);
                        elseif pairs(i,2) == currentVertex
                            newVertex = pairs(i,1);
                        end

                        if (newVertex ~= 0) && (newVertex ~= currentVertex)
                            if ismember(newVertex, coloredVertices)
                                % We already passed through this vertex, checking
                                % that coloring is consistent with previous choice
                                if colors(newVertex) ~= newColor
                                    coloring = [];
                                    return;
                                end
                            else
                                % We found a new vertex, let us color it
                                colors(newVertex) = newColor;
                                coColoredVertices = coColoredVertices + 1;
                                coloredVertices(coColoredVertices) = newVertex;
                            end
                        end
                    end
                end
            end
            ok = 1;
            coloring = [coloredVertices; colors(coloredVertices)];
        end
        
        function R = fromIndexMatrix(indexMatrix, generators, matrixType, field)
        % R = fromIndexMatrix(indexMatrix, generators, matrixType, field)
        % 
        % Creates an SDP matrix with additional structure. The produced
        % sdpvar matrix:
        %  - is invariant under the permutation group
        %  - satisfies the structure imposed by the index matrix: two
        %    matrix elements with same index are equal to each other
        % 
        % This is obtained by first performing an exact Reynolds
        % simplification of the matrix of indices so that it is invariant
        % under the group.
        %
        % The produced matrix is always invariant under hermitian
        % transposition.
        %
        % Args:
        %     indexMatrix: matrix with integer values corresponding to the
        %         label of the variable at each element. The actual index
        %         values are irrelevant.
        %     generators: a list of generators under which the matrix
        %         remains unchanged
        %     matrixType: one of the following:
        %         'full' : no particular structure
        %         'symmetric' : transpose-invariant matrix
        %         'hermitian' : hermitian matrix
        %     field: matrix elements are either 'real' or 'complex'
        %
        % Results:
        %     R: a CommutantVar object
        %
        % Example:
        %     indexMatrix = [1 1 3 4; 1 5 6 30; 3 6 10 11; 4 30 11 15];
        %     matrix = replab.CommutantVar.fromIndexMatrix(indexMatrix, {[4 1 2 3]}, 'symmetric', 'real')
        %
        % See also:
        %     replab.CommutantVar.fromPermutations
        %     replab.CommutantVar.fromSdpMatrix
        %     replab.CommutantVar.fromSymSdpMatrix

            % Basic tests
            assert(max(max(abs(indexMatrix - round(indexMatrix)))) == 0, 'The indexMatrix must be a matrix of integers.');
            d = size(indexMatrix, 1);
            for i = 1:length(generators)
                assert(length(generators{i}) == d, 'Generators and indexMatrix dimensions don''t match.');
            end
            assert(isequal(matrixType, 'full') || isequal(matrixType, 'symmetric') || isequal(matrixType, 'hermitian'), 'The matrix type must be ''full'', ''symmetric'' or ''hermitian''.');
            assert(isequal(field, 'real') || isequal(field, 'complex'), 'The field must be ''real'' or ''complex''.');
            if isequal(matrixType, 'hermitian') && isequal(field, 'real')
                matrixType = 'symmetric';
            end
            
            % First, we make the indexMatrix invariant under the group
            
            % We renumber all indices so they match default numbering
            [values, indices] = unique(indexMatrix(:), 'first');
            invPerm = sparse(values, 1, indices);
            indexMatrix = full(invPerm(indexMatrix));
            
            % We write the action of the generators on the matrix elements
            generators2 = cell(size(generators));
            M = reshape(1:d^2, [d d]);
            for i = 1:length(generators)
                generators2{i} = reshape(M(generators{i}, generators{i}), 1, d^2);
            end
            group = replab.Permutations(d^2).subgroup(generators2);

            % Identify the orbits
            orbits = group.orbits.blockIndex;
            
            % Make sure orbits are numbered in a similar way
            [values, indices] = unique(orbits(:), 'first');
            invPerm = sparse(values, 1, indices);
            orbits = full(invPerm(orbits));
            
            % Merge orbits with indices
            pairs = [reshape(indexMatrix, d^2, 1), orbits];

            if isequal(matrixType, 'symmetric')
                % Also impose that indexMatrix is symmetric
                pairs = [pairs; reshape(indexMatrix, d^2, 1) reshape(indexMatrix', d^2, 1)];
            end
            pairs = unique(sort(pairs, 2), 'rows');

            % We use a burning algorithm to identify all connected subsets
            subsets = replab.CommutantVar.burn(pairs);
            
            % We attribute the number to each element of each class
            images = zeros(max(unique(pairs)), 1);
            for i = 1:length(subsets)
                images(subsets{i}) = i;
            end
            
            % First we identify the index for each element
            % We substitute just one index per class of indices
            [values, indices] = unique(unique(indexMatrix(:)), 'first');
            invPerm = sparse(values, 1, indices);
            tmp = sparse(1:numel(indexMatrix), invPerm(reshape(indexMatrix,1,numel(indexMatrix))), true);
            %indexMatrix = reshape(tmp*images(values), d, d)
            
            % Now we can declare the desired number of variables and
            % attribute them
            if isequal(field, 'real')
                vars = sdpvar(length(subsets), 1);
            else
                % complex coefficients
                if isequal(matrixType, 'hermitian')
                    % find coefficients which must be real
                    imagesMatrix = reshape(tmp*images(values), d, d);
                    diagImages = unique(diag(imagesMatrix));
                    pairsH = [reshape(imagesMatrix, d^2, 1) reshape(imagesMatrix', d^2, 1)];
                    pairsH = unique(sort(pairsH, 2), 'rows');

                    subsetsH = replab.CommutantVar.burn(pairsH);
                    
                    % We distinguish between images which are real, and
                    % images which are conjugated to each other
                    subsetsR = {}; % real variables
                    subsetsC = {}; % mutually conjugated variables
                    for i = 1:length(subsetsH)
                        isReal = false;
                        for j = 1:length(diagImages)
                            if ismember(diagImages(j), subsetsH{i})
                                subsetsR{end+1} = subsetsH{i};
                                isReal = true;
                                break;
                            end
                        end
                        if ~isReal
                            % We try to find a 2-coloring of the graph
                            % corresponding to these indices. If it does
                            % not exist, then the variables must be real.
                            
                            % Fist we select the edges connecting vertices
                            % in this set
                            maskPairs = 0*pairsH;
                            for j = 1:length(subsetsH{i})
                                maskPairs = maskPairs + (pairsH == subsetsH{i}(j));
                            end
                            selPairs = pairsH(find(sum(maskPairs,2)),:);
                            
                            if sum(diff(selPairs,1,2) == 0) >= 1
                                % element must be equal to its conjugate,
                                % so it must be real
                                subsetsR{end+1} = subsetsH{i};
                            else
                                % We try to split the element with respect
                                % to conjugacy
                                [colorable, coloring] = replab.CommutantVar.graphIsBipartite(selPairs);
                                if colorable == 0
                                    subsetsR{end+1} = subsetsH{i};
                                else
                                    subsetsC{end+1} = coloring;
                                end
                            end
                        end
                    end
                    
                    % the real variables
                    nbRealVariables = length(subsetsR);
                    varsR = sdpvar(nbRealVariables,1);
                    
                    % placement of the real variables
                    indR1 = [subsetsR{:}];
                    indR2 = zeros(size(indR1));
                    co = 1;
                    for i = 1:nbRealVariables
                        indR2(co:co+length(subsetsR{i})-1) = i;
                        co = co + length(subsetsR{i});
                    end
                    
                    % the complex variables
                    nbComplexVariables = length(subsetsC);
                    varsC = sdpvar(nbComplexVariables,1,'full','complex');
                    varsC = reshape([varsC conj(varsC)].', 2*nbComplexVariables,1);
                    
                    % placement of the complex variables
                    indC1 = [subsetsC{:}];
                    indC2 = zeros(1,size(indC1,2));
                    co = 1;
                    for i = 1:nbComplexVariables
                        indC2(co+find(subsetsC{i}(2,:)==1)-1) = 2*(i-1)+1;
                        indC2(co+find(subsetsC{i}(2,:)==2)-1) = 2*(i-1)+2;
                        co = co + size(subsetsC{i},2);
                    end
                    
                    % All together
                    if nbComplexVariables == 0
                        ind1 = indR1;
                        ind2 = indR2;
                        field = 'real';
                        if isequal(matrixType, 'hermitian')
                            matrixType = 'symmetric';
                        end
                    else
                        ind1 = [indR1, indC1(1,:)];
                        ind2 = [indR2, indR2(end)+indC2];
                    end
                    
                    vars = sparse(ind1, ind2, 1)*[varsR; varsC];
                else
                    vars = sdpvar(length(subsets), 1, 'full', 'complex');
                end
            end
            sdpMatrix = reshape(tmp*vars(images(values)), d, d);
            
            R = replab.CommutantVar(generators, sdpMatrix, 1, matrixType, field);
        end

        function R = fromSdpMatrix(sdpMatrix, generators)
        % R = fromSdpMatrix(sdpMatrix, generators)
        % 
        % Imposes symmetry constraints onto an existing SDP matrix. This
        % creates a CommutantVar matrix which satisfies both:
        %     - the structure encoded into sdpMatrix
        %     - invariance under joint permutations of its lines and
        %       columns by the provided generators.
        % The type of matrix (full/symmetric/hermitian) as well as the
        % field (real/complex) is inferred from the provided matrix.
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
        %     replab.CommutantVar.fromIndexMatrix
        %     replab.CommutantVar.fromSymSdpMatrix

            if issymmetric(sdpMatrix)
                matrixType = 'symmetric';
            elseif ishermitian(sdpMatrix)
                matrixType = 'hermitian';
            else
                matrixType = 'full';
            end
            
            if isreal(sdpMatrix)
                field = 'real';
            else
                field = 'complex';
            end
            
            R = replab.CommutantVar(generators, sdpMatrix, 0, matrixType, field);
        end

        function R = fromSymSdpMatrix(sdpMatrix, generators)
        % R = fromSymSdpMatrix(sdpMatrix, generators)
        % 
        % Block-diagonalizes an existing SDP matrix which is already
        % invariant under the group generators.
        % The type of matrix (full/symmetric/hermitian) as well as the
        % field (real/complex) is inferred from the provided matrix.
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
        %     replab.CommutantVar.fromIndexMatrix
        %     replab.CommutantVar.fromSdpMatrix

            if issymmetric(sdpMatrix)
                matrixType = 'symmetric';
            elseif ishermitian(sdpMatrix)
                matrixType = 'hermitian';
            else
                matrixType = 'full';
            end
            
            if isreal(sdpMatrix)
                field = 'real';
            else
                field = 'complex';
            end
            
            R = replab.CommutantVar(generators, sdpMatrix, 1, matrixType, field);
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
        
        % TODO: integrate this function into the interaction with the Str
        % class
        
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
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]}, 'symmetric, 'real')
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
        %     matrix = replab.CommutantVar.fromPermutations({[3 1 2]}, 'symmetric, 'real')
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
        %     M: matrix
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     matrix.blockMask
        
            M = zeros(self.dim, self.dim);
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

        function M = fullBlockMatrix(self)
        % M = fullBlockMatrix(self)
        %
        % Returns the full matrix in its block-diagonal form.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     M: sdpvar matrix in block diagonal form
        %
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     see(matrix.fullBlockMatrix)
            
            if isempty(self.fullBlockMatrix_)
                % We construct the matrix for the first time
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
                self.fullBlockMatrix_ = M;
            else
                M = self.fullBlockMatrix_;
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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     see(matrix.fullMatrix)
        
            if ~isempty(self.sdpMatrix_)
                M = self.sdpMatrix_;
            else
                instructions = struct('type', '()');
                instructions.subs = {':', ':'};
                M = subsref(self, instructions);
            end
        end

        function M = U(self)
        % M = U(self)
        %
        % Returns the block-diagonalizing unitary.
        %
        % Args:
        %     self: CommutantVar object
        %
        % Returns:
        %     U: matrix
        % 
        % Example:
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     matrix.U
        
            M = full(self.U_);
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
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
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
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
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

        function lastIndex = end(self, k, n)
        % end - returns the last index
        %
        % lastIndex = replab.CommutantVar.end(self, k, n)
        % Returns the last index in an indexing expression such as 
        % self(1,2:end) or self(end).
        %
        % Args:
        %     self: CommutantVar object
        %     k: 
        %     n: number of dimensions on which the indexing is done
        %
        % Returns:
        %     lastIndex: double
        %
        % See also:
        %     replab.CommutantVar.subsref

            if n == 1
                % Then the indexing is done with only one index, as in a(end)
                lastIndex = prod(size(self));
            else
                % Then the indexing is done with two indices, as in a(end,1:2)
                if (k < 1) || (k > 2)
                    error('Only two dimensions available');
                end
                s = size(self);
                lastIndex = s(k);
            end
        end

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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     matrix(1)
        %
        % See also:
        %     subsref
        
            switch varargin{1}(1).type
                case '()'
                    % If the matrix in full form is known, we extract the
                    % elements directly from it
                    if ~isempty(self.sdpMatrix_)
                        [varargout{1:nargout}] = subsref(self.sdpMatrix_, varargin{1});
                        return;
                    end
                    
                    % The matrix in full form has not been constructed, we
                    % rely on the matrix in block form
                    switch length(varargin{1}.subs)
                        case 1
                            % call of the form self([1 2 4 5])
                            [varargout{1:nargout}] = subsref(self.fullMatrix, varargin{1});
                        case 2
                            % call of the form self([1 2], [4 5])
                            
                            % If needed, we replace ':' by actual indices
                            for i = 1:2
                                if ischar(varargin{1}.subs{i}) && isequal(varargin{1}.subs{i}, ':')
                                    varargin{1}.subs{i} = 1:size(self,i);
                                end
                            end
                            
                            % Now we extract only the requested part
                            M = self.fullBlockMatrix;
                            if length(varargin{1}.subs{1})*size(self,2) <= size(self,1)*length(varargin{1}.subs{2})
                                varargout{1} = (self.U_(varargin{1}.subs{1},:)*M)*(self.U_(varargin{1}.subs{2},:)');
                            else
                                varargout{1} = self.U_(varargin{1}.subs{1},:)*(M*(self.U_(varargin{1}.subs{2},:)'));
                            end
                            
                            % Make sure symmetry/hermiticity is preserved
                            if isequal(varargin{1}.subs{1}, varargin{1}.subs{2})
                                if isequal(self.matrixType, 'symmetric')
                                    varargout{1} = (varargout{1} + varargout{1}.')/2;
                                elseif isequal(self.matrixType, 'hermitian')
                                    varargout{1} = (varargout{1} + varargout{1}')/2;
                                end
                            end
                            
                        otherwise
                            error('Too many indices for 2-dimensional object');
                    end
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
        %     matrix1 = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     matrix2 = replab.CommutantVar.fromPermutations({[2 1 3]}, 'symmetric, 'real')
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
                    rotatedY = X.U_'*Y.fullMatrix*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
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
                rotatedY = X.U_'*Y*X.U_;

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
                    rotatedY = X.U_'*Y.fullMatrix*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
                end
    
                % We want to make sure that the symmetry or hermitianity is
                % preserved if possible. So we check what symmetry is
                % expected in the result of the operation
                if isa(Y, 'replab.CommutantVar')
                    YType = Y.matrixType;
                    YField = Y.field;
                    
                    YSdpMatrix = Y.sdpMatrix_;
                else
                    if isreal(Y)
                        YField = 'real';
                    else
                        YField = 'complex';
                    end
                    if issymmetric(Y)
                        YType = 'symmetric';
                    elseif ishermitian(Y)
                        YType = 'hermitian';
                    else
                        YType = 'full';
                    end
                    
                    YSdpMatrix = Y;
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
                    
                    % We make sure that the symmetry or hermitianity is
                    % preserved if possible
                    if ~isequal(YType, 'full')
                        if isequal(YType, 'symmetric')
                            shouldBeZero = constantBlock - constantBlock.';
                        elseif isequal(YType, 'hermitian')
                            shouldBeZero = constantBlock - constantBlock';
                        end
                        if isa(shouldBeZero, 'sdpvar')
                            indices = [0 getvariables(shouldBeZero)];
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Non-', YType,' component of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), ' introduced.']);
                                end
                                maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Non-', YType,' component of ', num2str(max(max(abs(getbasematrix(shouldBeZero))))), ' introduced.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                        end
                        if isequal(YType, 'symmetric') && ~issymmetric(constantBlock)
                            % We force exact symmetry
                            constantBlock = (constantBlock + constantBlock.')/2;
                        elseif isequal(YType, 'hermitian') && ~ishermitian(constantBlock)
                            % We force exact hermiticity
                            constantBlock = (constantBlock + constantBlock')/2;
                        end
                    end
                    Z.blocks{i} = Z.blocks{i} + constantBlock;
                    co = co + d*size(X.blocks{i}, 1);
                end
                
                % We keep track of the linear constraints
                if ~isempty(X.sdpMatrix_) && ~isempty(YSdpMatrix)
                    Z.sdpMatrix_ = X.sdpMatrix_ + YSdpMatrix;
                else
                    Z.sdpMatrix_ = [];
                end
                if isa(Y, 'replab.CommutantVar')
                    Z.linearConstraints = [X.linearConstraints, Y.linearConstraints];
                else
                    Z.linearConstraints = X.linearConstraints;
                end
                
                % We update the matrix attributes
                newType = 'full';
                if isequal(X.matrixType, 'symmetric') && isequal(YType, 'symmetric')
                    newType = 'symmetric';
                elseif isequal(X.matrixType, 'hermitian') && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif (isequal(X.field, 'real') && isequal(X.matrixType, 'symmetric')) && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif isequal(X.matrixType, 'hermitian') && (isequal(YField, 'real') && isequal(YType, 'symmetric'))
                    newType = 'hermitian';
                end
                Z.matrixType = newType;
                if isequal(X.field, 'real') && isequal(YField, 'real')
                    Z.field = 'real';
                else
                    Z.field = 'complex';
                end
                
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg + CommutantVar
                Z = Y+X;
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end

            
            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if isequal(newType, 'symmetric') && (maxNonHermiticity > epsilonWarning)
                warning(['Non-symmetry of order ', num2str(maxEpsilonFound), ' was corrected.']);
            elseif isequal(newType, 'hermitian') && (maxNonHermiticity > epsilonWarning)
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
                    rotatedY = X.U_'*Y.fullMatrix*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
                end
    
                % We want to make sure that the symmetry or hermitianity is
                % preserved if possible. So we check what symmetry is
                % expected in the result of the operation
                if isa(Y, 'replab.CommutantVar')
                    YType = Y.matrixType;
                    YField = Y.field;
                    
                    YSdpMatrix = Y.sdpMatrix_;
                else
                    if isreal(Y)
                        YField = 'real';
                    else
                        YField = 'complex';
                    end
                    if issymmetric(Y)
                        YType = 'symmetric';
                    elseif ishermitian(Y)
                        YType = 'hermitian';
                    else
                        YType = 'full';
                    end
                    
                    YSdpMatrix = Y;
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
                    
                    % We make sure that the symmetry or hermitianity is
                    % preserved if possible
                    if ~isequal(YType, 'full')
                        if isequal(YType, 'symmetric')
                            shouldBeZero = constantBlock - constantBlock.';
                        elseif isequal(YType, 'hermitian')
                            shouldBeZero = constantBlock - constantBlock';
                        end
                        if isa(shouldBeZero, 'sdpvar')
                            indices = [0 getvariables(shouldBeZero)];
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Non-', YType,' component of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), ' introduced.']);
                                end
                                maxNonHermiticity = max([maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind))))]);
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Non-', YType,' component of ', num2str(max(max(abs(getbasematrix(shouldBeZero))))), ' introduced.']);
                            end
                            maxNonHermiticity = max([maxNonHermiticity, max(max(abs(shouldBeZero)))]);
                        end
                        if isequal(YType, 'symmetric') && ~issymmetric(constantBlock)
                            % We force exact symmetry
                            constantBlock = (constantBlock + constantBlock.')/2;
                        elseif isequal(YType, 'hermitian') && ~ishermitian(constantBlock)
                            % We force exact hermiticity
                            constantBlock = (constantBlock + constantBlock')/2;
                        end
                    end
                    Z.blocks{i} = Z.blocks{i} - constantBlock;
                    co = co + d*size(X.blocks{i}, 1);
                end
                
                % We keep track of the linear constraints
                if ~isempty(X.sdpMatrix_) && ~isempty(YSdpMatrix)
                    Z.sdpMatrix_ = X.sdpMatrix_ - YSdpMatrix;
                else
                    Z.sdpMatrix_ = [];
                end
                if isa(Y, 'replab.CommutantVar')
                    Z.linearConstraints = [X.linearConstraints, Y.linearConstraints];
                else
                    Z.linearConstraints = X.linearConstraints;
                end
                
                % We update the matrix attributes
                newType = 'full';
                if isequal(X.matrixType, 'symmetric') && isequal(YType, 'symmetric')
                    newType = 'symmetric';
                elseif isequal(X.matrixType, 'hermitian') && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif (isequal(X.field, 'real') && isequal(X.matrixType, 'symmetric')) && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif isequal(X.matrixType, 'hermitian') && (isequal(YField, 'real') && isequal(YType, 'symmetric'))
                    newType = 'hermitian';
                end
                Z.matrixType = newType;
                if isequal(X.field, 'real') && isequal(YField, 'real')
                    Z.field = 'real';
                else
                    Z.field = 'complex';
                end
            elseif ~isa(X, 'replab.CommutantVar') && isa(Y, 'replab.CommutantVar')
                % sthg - CommutantVar
                Z = -(Y-X);
            else
                error('Neither of the two arguments is of type replab.CommutantVar');
            end


            % We produce a warning if some small but not too small coefficients
            % have been neglected
            if isequal(newType, 'symmetric') && (maxNonHermiticity > epsilonWarning)
                warning(['Non-symmetry of order ', num2str(maxEpsilonFound), ' was corrected.']);
            elseif isequal(newType, 'hermitian') && (maxNonHermiticity > epsilonWarning)
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
            X.sdpMatrix_ = -X.sdpMatrix_;
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
                
                % Keeping track of linear constraints
                if ~isempty(X.sdpMatrix_)
                    Z.sdpMatrix_ = X.sdpMatrix_.*Y;
                end
                
                % Extract attributes of Y
                if isreal(Y)
                    YField = 'real';
                else
                    YField = 'complex';
                end
                if issymmetric(Y)
                    YType = 'symmetric';
                elseif ishermitian(Y)
                    YType = 'hermitian';
                else
                    YType = 'full';
                end
                
                % We update the matrix attributes
                newType = 'full';
                if isequal(X.matrixType, 'symmetric') && isequal(YType, 'symmetric')
                    newType = 'symmetric';
                elseif isequal(X.matrixType, 'hermitian') && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif (isequal(X.field, 'real') && isequal(X.matrixType, 'symmetric')) && isequal(YType, 'hermitian')
                    newType = 'hermitian';
                elseif isequal(X.matrixType, 'hermitian') && (isequal(YField, 'real') && isequal(YType, 'symmetric'))
                    newType = 'hermitian';
                end
                Z.matrixType = newType;
                if isequal(X.field, 'real') && isequal(YField, 'real')
                    Z.field = 'real';
                else
                    Z.field = 'complex';
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
            
            % Keeping track of linear constraints
            if ~isempty(X.sdpMatrix_)
                Z.sdpMatrix_ = X.sdpMatrix_./Y;
            end

            % Extract attributes of Y
            if isreal(Y)
                YField = 'real';
            else
                YField = 'complex';
            end

            % We update the matrix attributes
            newType = 'full';
            if isequal(X.matrixType, 'symmetric')
                newType = 'symmetric';
            elseif isequal(X.matrixType, 'hermitian') && isequal(YField, 'real')
                newType = 'hermitian';
            end
            Z.matrixType = newType;
            if isequal(X.field, 'real') && isequal(YField, 'real')
                Z.field = 'real';
            else
                Z.field = 'complex';
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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
        %     trace(matrix)
        % 
        % See also:
        %     trace

            if ~isempty(self.sdpMatrix_)
                X = diag(self.sdpMatrix_);
            else
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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
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
                    rotatedY = X.U_'*Y.U_*Y.blockMask*Y.U_'*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
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
                    rotatedY = X.U_'*Y.U_*Y.blockMask*Y.U_'*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
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
        %     matrix = replab.CommutantVar.fromPermutations({[2 3 1]}, 'symmetric, 'real')
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
                    rotatedY = X.U_'*Y.U_*Y.blockMask*Y.U_'*X.U_;
                else
                    rotatedY = X.U_'*Y*X.U_;
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
