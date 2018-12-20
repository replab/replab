classdef SymmetricSDP
% A block-diagonalizing preprocessor for SDPs in SeDuMi format.
%
% The problem data is stored in (A, b, c, K) according to the SeDuMi
% documentation, with size(A, 2) the number of dual variables.
%
% We only support a single SDP block, and no other cones
% (so K.f = K.l = 0, K.q = K.r = [] and K.s = s).
%
% We provide the symmetry group acting on the SDP block, such that
% rho(g) * X * rho(g)' = X for all elements g in the group G.
%
% The finite group G and its representation rho are jointly written as follows.
%
% Without loss of generality, a finite group can be represented by permutations.
% Then G is a (row) cell array containing the N generators of this permutation
% group.
    
% rho is also a 1xnG cell array, with each element describing a
% s x s double matrix, such that rho{i} is the image of the generator G{i}.

    properties
        A; % sedumi data
        b; % sedumi data
        c; % sedumi data
        s; % size of single SDP block present
        m; % number of dual variables
        K; % sedumi data
        G; % cell array of generators as permutations
        rho; % cell array of generator images defining the representation
        rep; % computed replab.FiniteGroupRep
    end
    
    methods
        
        function self = SymmetricSDP(A, b, c, K, G, rho)
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
            group = replab.PermutationGroup.fromGenerators(G{:});
            % check for unitarity
            tol = 1e-12;
            for i = 1:length(rho)
                assert(norm(rho{i}*rho{i}' - eye(self.s)) < tol);
            end
            self.rep = group.representation(self.s, rho, true);
        end
        
        function vec1 = project(self, vec)
            rep = self.rep;
            s = self.s;
            nc = s*s;
            nb1 = rep.irreducible.nComponents; % number of blocks
            s1 = rep.irreducible.centralizerBlockDimensions; % output block sizes
            nc1 = sum(s1.^2);
            vec1 = zeros(nc1, 1);
            shift = 0;
            M = reshape(vec, [s s]);
            for r = 1:nb1 % iterate over represenattions
                % apply Reynolds
                block = rep.irreducible.blockOfCentralizer(M, r, 0); 
                % store block flattened 
                nels = prod(size(block));
                vec1(shift+(1:nels)) = block(:);
                shift = shift + nels;
            end
        end
        
        function [A1 b1 c1 K1] = blockDiagonalize(self)
            rep = self.rep;
            m = self.m; % number of dual variables
            s1 = rep.irreducible.centralizerBlockDimensions; % output block sizes
            nc1 = sum(s1.^2);
            b1 = self.b;
            c1 = self.project(self.c);
            A1 = zeros(nc1, m);
            for i = 1:m
                A1(:,i) = self.project(self.A(:,i));
            end
            K1 = struct('f', 0, 'l', 0, 'q', 0, 'r', 0, 's', s1);
            newData = struct('K', K1, 'A', A1, 'b', self.b, 'c', c1);            
        end
        
    end
end


