classdef SedumiData
% A block-diagonalizing preprocessor for SDPs in SeDuMi format.
%
% The problem data is stored in (A, b, c, K) according to the SeDuMi
% documentation, with size(A, 2) the number of dual variables.
%
% We only support a single SDP block, and no other cones
% (so K.f = K.l = 0, K.q = K.r = [] and K.s = s).
%
% We provide the symmetry group acting on the SDP block X, such that
% rho(g) * X * rho(g)' = X for all elements g in the group G.
%
% The finite group G and its representation rho are jointly written as follows.
%
% Without loss of generality, a finite group can be represented by permutations.
% Then G is a (row) cell array containing the N generators of this permutation
% group.
    
% rho is also a 1xnG cell array, with each element describing a
% s x s double matrix, such that rho{i} is the image of the generator G{i}.
%
% The representation must be unitary.

    properties
        A % double: Data matrix of size n x m, where n is the number of primal scalar variables, 
          %         and m the nubmer of primal constraints (Sedumi data)
        b % double: Constraint right hand side (Sedumi data)
        c % double: Primal objective (Sedumi data)
        K % struct: Cone specification (Sedumi data)
        m % number of dual variables
        s % size of single SDP block present
        G % cell array of generators as permutations
        rho % cell array of generator images defining the representation
        rep % replab.Rep: group representation commuting with the SDP block
    end
    
    methods
        
        function self = SedumiData(A, b, c, K, G, rho)
            assert(~isempty(G), 'Can only process SDPs with symmetries, but no generators present');
            self.A = A;
            self.b = b;
            self.c = c;
            self.K = K;
            assert(K.f == 0, 'No support for free cone');
            assert(K.l == 0, 'No support for nonneg cone');
            assert(length(K.q) == 0, 'No support for quadratic cones');
            assert(length(K.r) == 0, 'No support for quadratic cones');
            assert(length(K.s) == 1, 'Can only process a single SDP block');
            self.s = K.s;
            self.m = size(A, 2);
            self.G = G;
            self.rho = rho;
            % build permutation group
            assert(~isempty(G), 'Must be nontrivial permutation group');
            ds = length(G{1}); % group domain size
            group = replab.Permutations(ds).subgroup(G);
            % check for unitarity
            tol = 1e-12;
            for i = 1:length(rho)
                assert(norm(rho{i}*rho{i}' - eye(self.s)) < tol);
            end
            self.rep = group.repByImages('R', self.s, true, rho);
        end
        
        function vec1 = project(self, vec)
            rep = self.rep;
            I = rep.decomposition;
            s = self.s;
            nc = s*s;
            nb1 = I.nComponents; % number of blocks
            C = I.asRep.commutant;
            s1 = arrayfun(@(c) C.singleBlockSize(c), 1:I.nComponents);
            nc1 = sum(s1.^2);
            vec1 = zeros(nc1, 1);
            shift = 0;
            M = reshape(vec, [s s]);
            for r = 1:nb1 % iterate over representations
                          % apply Reynolds
                block = C.blockOfParentElement(M, r, false);
                % store block flattened 
                nels = prod(size(block));
                vec1(shift+(1:nels)) = block(:);
                shift = shift + nels;
            end
            assert(length(vec1) == nc1);
        end
        
        function [A1 b1 c1 K1] = blockDiagonalize(self)
            rep = self.rep;
            I = rep.decomposition;
            m = self.m; % number of dual variables
            C = I.asRep.commutant;
            s1 = arrayfun(@(c) C.singleBlockSize(c), 1:I.nComponents); % output block sizes
            nc1 = sum(s1.^2);
            b1 = self.b;
            c1 = self.project(self.c);
            A1 = zeros(nc1, m);
            for i = 1:m
                p = self.project(self.A(:,i));
                A1(:,i) = p;
            end
            K1 = struct('f', 0, 'l', 0, 'q', 0, 'r', 0, 's', s1);
        end
        
        function saveBlockDiagonalizedMatFile(self, filename)
            [A1 b1 c1 K1] = self.blockDiagonalize;
            newData = struct('K', K1, 'A', A1, 'b', b1, 'c', c1);
            save(filename, '-struct', 'newData');
        end
        
    end
    
    methods (Static)
        
        function sedumiData = fromMatFile(filename)
            data = load(filename);
            if ~isfield(data, 'G') || ~isfield(data, 'rho')
                error('The file does not contain symmetry information');
            end
            sedumiData = replab.SedumiData(data.A', data.b, data.c, data.K, data.G, data.rho);
        end
        
    end
    
end
