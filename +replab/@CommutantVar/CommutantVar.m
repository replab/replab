classdef CommutantVar < replab.Str
% A matrix variable satisfying some symmetry constraints
%
% example: replab.CommutantVar.fromPermutations({[3 1 2]})

% Current limitations:
% - The sdp matrix produced is currently always square, real, and symmetric
% - Only supports complex and quaternionic representations with dimension 2
%   and 4 respectively.
% - limited arithmetic, essentially we can just:
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
        
        function self = CommutantVar(U, dimensions1, multiplicities, types)
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
        
        function R = fromPermutations(generators)
            assert(iscell(generators), 'Please specify generators in cell array.');
            n = size(generators{1}, 2);
            group = replab.SignedPermutations(n).subgroup(generators);
            
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
                        
            R = replab.CommutantVar(U, dimensions1, multiplicities, types);
        end
        
    end

end
