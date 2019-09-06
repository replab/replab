classdef CommutantVar < replab.Str
%classdef (InferiorClasses = {?sdpvar,?gem,?sgem}) CommutantVar < replab.Str %& matlab.mixin.Copyable % Use this declaration for correct class precedence in matlab
%
% A matrix variable satisfying some symmetry constraints
%
% example: replab.CommutantVar.fromPermutations({[3 1 2]})

% Warning: this object inherits from a handle object, therefore it is also
% a handle object. To copy this object use the 'copy' method to obtain two
% identical but independent objects.

% Current limitations:
% - The sdp matrix produced is currently always square, real, and symmetric
% - Only supports complex and quaternionic representations with dimension 2
%   and 4 respectively.

    properties (SetAccess = protected)
        U; % Unitary operator block-diagonalizing the matrix
        nComponents; % The block structure : dimension1 x multiplicity
        dimensions1;
        multiplicities;
        types; % The representation type
        dim; % matrix dimension
        blocks; % The sdp blocks corresponding to each irreducible representation
    end

    properties (SetAccess = public)
        sdpMatrix; % An sdpMatrix that the CommutantVar must match
    end

    methods


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                          inherited classes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function s = headerStr(self)
            s = sprintf('Commutant variable %dx%d (%d blocks, %d scalar variables)', self.dim, self.dim, length(self.blocks), self.nbVars());
        end

        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'blocks';
            names{1, end+1} = 'nComponents';
            names{1, end+1} = 'sdpMatrix';
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.Str(self);

            names{1, end+1} = 'blocks';
            dims = zeros(1,length(self.blocks));
            for i = 1:length(self.blocks)
                dims(i) = size(self.blocks{i},1);
            end
            values{1, end+1} = dims;

            if ~isempty(self.sdpMatrix)
                names{1, end+1} = 'sdpMatrix';
                values{1, end+1} = self.sdpMatrix;
            end
        end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                               constructors
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function self = CommutantVar(U, dimensions1, multiplicities, types, sdpMatrix)
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

            % We keep a copy of the SDP matrix constraint if provided
            self.sdpMatrix = sdpMatrix;
        end

        function R = copy(rhs)
            % copy object

            % Create a new object
            R = replab.CommutantVar(rhs.U, rhs.dimensions1, rhs.multiplicities, rhs.types, rhs.sdpMatrix);

            % Override properties
%             fns = properties(rhs);
%             for i=1:length(fns)
%                 R.(fns{i}) = rhs.(fns{i});
%             end
            % The above is not supported on octave, so we copy all elements
            % by hand...
            R.U = rhs.U;
            R.nComponents = rhs.nComponents;
            R.dimensions1 = rhs.dimensions1;
            R.multiplicities = rhs.multiplicities;
            R.types = rhs.types;
            R.dim = rhs.dim;
            R.blocks = rhs.blocks;
            R.sdpMatrix = rhs.sdpMatrix;
        end

    end

    methods (Static) % Factory methods

        function R = fromPermutations(generators)
            R = replab.CommutantVar.fromSDPMatrix([], generators);
        end

        function R = fromSDPMatrix(sdpMatrix, generators)
            assert(nargin >= 2, 'Not enough arguments.');

            assert(iscell(generators), 'Please specify generators in cell array.');
            n = size(generators{1}, 2);
            group = replab.SignedPermutations(n).subgroup(generators);

            if ~isempty(sdpMatrix)
                assert(isequal(size(sdpMatrix), [n n]), 'Wrong matrix or group dimension.');
            end

            irrDecomp = group.naturalRep.decomposition;
            U = zeros(n, 0);
            dimensions1 = [];
            multiplicities = [];
            types = '';
            for i = 1:irrDecomp.nComponents
                component = irrDecomp.component(i);
                dimensions1 = [dimensions1 component.copyDimension];
                multiplicities = [multiplicities component.multiplicity];
                types(i) = component.copy(1).realDivisionAlgebra.shortName;
                for j = 1:component.multiplicity
                    copy = component.copy(j);
                    % correction, as the new RepLAB convention
                    % is to store basis vectors as row vectors
                    U = [U copy.U'];
                end
            end

            R = replab.CommutantVar(U, dimensions1, multiplicities, types, sdpMatrix);
        end

    end


    methods

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                       access to class properties
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

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

        function [s1 s2] = size(self, d)
            % returns the matrix size
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
        % Returns the desired block
            M = self.blocks{i};
        end

        function M = blockMask(self)
        % Returns a mask showing the block structure
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
        % Constructs the full SDP matrix in the natural basis
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
        % Returns the SDP variable used be the object
            vars = getvariables(self.blocks{1});
            for i = 2:self.nComponents
                vars = [vars, getvariables(self.blocks{i})];
            end
            vars = unique(vars);
        end

        function n = nbVars(self)
        % Returns the number of SDP variable used be the object
            n = 0;
            for i = 1:self.nComponents
                n = n + length(getvariables(self.blocks{i}));
            end
        end

        function val = value(self)
        % Returns the current numerical value of the object
            val = value(self.fullMatrix);
        end

        function see(self)
        % Displays internal info about the matrix composition. For now, we
        % delegate to yalmip
            see(self.fullMatrix);
        end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                            utility methods
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

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

        function okLevel = compatibleWith(X, Y)
        % checks whether the block structure of Y is compatible with that of X
        %
        % the output is:
        %  0 if Y is not compatible with X
        %  1 if Y has the block structure of X
        %  2 if Y has the block structure of X and identical blocks in X
        %    correspond to identical blocks in Y

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track if some rounding is done
            maxOuterEpsilonFound = 0;
            maxEpsilonFound = 0;
            epsilonWarning = 1e-14;

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
                    indices = getvariables(shouldBeZero);
                    for ind = indices
                        if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                            okLevel = 0;
                            break;
                        end
                        maxOuterEpsilonFound = max(maxOuterEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                    end

                    % Now we check the intro-block structure
                    if okLevel == 2
                        % Check full compatibility
                        ideal = rotatedY(co + (1:d:d*size(X.blocks{i},1)), co + (1:d:d*size(X.blocks{i},2)));
                        ideal = kron(ideal, eye(d));
                        shouldBeZero = ideal - rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2)));

                        % we check if the coefficients are negligible
                        indices = getvariables(shouldBeZero);
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                okLevel = 1;
                                maxEpsilonFound = 0;
                                break;
                            end
                            maxEpsilonFound = max(maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                        end
                    end
                    if okLevel == 1
                        % Check partial compatibility
                        mask = kron(ones(size(X.blocks{i})), eye(d));
                        shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2))).*(1-mask);

                        % we check if the coefficients are negligible
                        indices = getvariables(shouldBeZero);
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                okLevel = 0;
                                return;
                            end
                            maxEpsilonFound = max(maxEpsilonFound, max(max(abs(getbasematrix(shouldBeZero,ind)))));
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
                    maxOuterEpsilonFound = max(maxOuterEpsilonFound, max(max(abs(shouldBeZero))));

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
                        maxEpsilonFound = max(maxEpsilonFound, max(max(abs(shouldBeZero))));
                    end
                    if okLevel == 1
                        % Check partial compatibility
                        mask = kron(ones(size(X.blocks{i})), eye(d));
                        shouldBeZero = rotatedY(co + (1:d*size(X.blocks{i},1)), co + (1:d*size(X.blocks{i},2))).*(1-mask);
                        if max(max(abs(shouldBeZero))) > epsilon
                            okLevel = 0;
                            return;
                        end
                        maxEpsilonFound = max(maxEpsilonFound, max(max(abs(shouldBeZero))));
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
            if max(maxOuterEpsilonFound, maxEpsilonFound) > epsilonWarning
                warning(['Block structure mismatch by ', num2str(maxEpsilonFound), ' ignored.']);
            end
        end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                             algebra methods
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        function Z = plus(X,Y)
        % addition operator

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-14;

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
                        indices = getvariables(shouldBeZero);
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                            end
                            maxNonHermiticity = max(maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                        end
                    else
                        if max(max(abs(shouldBeZero))) > epsilon
                            error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                        end
                        maxNonHermiticity = max(maxNonHermiticity, max(max(abs(shouldBeZero))));
                    end
                    if ~ishermitian(constantBlock)
                        % We force exact hermiticity
                        constantBlock = (constantBlock + constantBlock')/2;
                    end
                    Z.blocks{i} = Z.blocks{i} + constantBlock;
                    co = co + d*size(X.blocks{i},1);
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
        % substraction operator

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-14;

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
                        indices = getvariables(shouldBeZero);
                        for ind = indices
                            if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                            end
                            maxNonHermiticity = max(maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                        end
                    else
                        if max(max(abs(shouldBeZero))) > epsilon
                            error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                        end
                        maxNonHermiticity = max(maxNonHermiticity, max(max(abs(shouldBeZero))));
                    end
                    if ~ishermitian(constantBlock)
                        % We force exact hermiticity
                        constantBlock = (constantBlock + constantBlock')/2;
                    end
                    Z.blocks{i} = Z.blocks{i} - constantBlock;
                    co = co + d*size(X.blocks{i},1);
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
            % unary minus operator

            X = self.copy;
            for i = 1:X.nComponents
                X.blocks{i} = -X.blocks{i};
            end
        end

        function Z = times(X, Y)
            % element-wise multiplication

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
            % element-wise right division

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
            % element-wise left division

            Z = Y./X;
        end

        function Z = mtimes(X, Y)
            % matrix multiplication

            % Check that dimensions are compatible
            if (numel(X) ~= 1) && (numel(Y) ~= 1)
                error('Incompatible size for multiplication, use fullMatrix for non-scalar multiplications.');
            end

            Z = X.*Y;
        end

        function Z = mrdivide(X, Y)
            % matrix right division

            % Only scalar division is supported
            if numel(Y) ~= 1
                error('Use fullMatrix for non-scalar division.');
            end

            Z = X./Y;
        end

        function Z = mldivide(X, Y)
            % matrix left division

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
        % trace operator
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

        function F = ge(X,Y)
            % greater or equal constraint

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-14;

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar') && (numel(Y) == 1)
                % CommutantVar >= scalar

                F = (X.blocks{1} >= Y);
                for i = 2:X.nComponents
                    F = [F, X.blocks{i} >= Y];
                end
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
                            indices = getvariables(shouldBeZero);
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                                end
                                maxNonHermiticity = max(maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                            end
                            maxNonHermiticity = max(maxNonHermiticity, max(max(abs(shouldBeZero))));
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
            % greater or equal constraint

            % Numerical tolerance to decide whether numbers are close to zero in
            % this function.
            epsilon = 1e-10;

            % We keep track of the encountered non-hermiticities
            maxNonHermiticity = 0;
            epsilonWarning = 1e-14;

            % We examine each case independently
            if isa(X, 'replab.CommutantVar') && ~isa(Y, 'replab.CommutantVar') && (numel(Y) == 1)
                % CommutantVar <= scalar

                F = (X.blocks{1} <= Y);
                for i = 2:X.nComponents
                    F = [F, X.blocks{i} <= Y];
                end
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
                            indices = getvariables(shouldBeZero);
                            for ind = indices
                                if max(max(abs(getbasematrix(shouldBeZero,ind)))) > epsilon
                                    error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(getbasematrix(shouldBeZero,ind))))), '.']);
                                end
                                maxNonHermiticity = max(maxNonHermiticity, max(max(abs(getbasematrix(shouldBeZero,ind)))));
                            end
                        else
                            if max(max(abs(shouldBeZero))) > epsilon
                                error(['Positivity requires hermitian matrices. Non-hermiticity of ', num2str(max(max(abs(shouldBeZero)))), '.']);
                            end
                            maxNonHermiticity = max(maxNonHermiticity, max(max(abs(shouldBeZero))));
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
            % greater than constraint
            warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
            F = (X >= Y);
        end

        function F = lt(X,Y)
            % lesser than constraint
            warning('Strict inequalities cannot be guaranteed. Imposing a non-strict constraint instead.');
            F = (X <= Y);
        end

    end

end
