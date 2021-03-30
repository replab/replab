run ~/w/replab/replab_init.m
addpath(genpath('~/software/QETLAB'));
% We compute the inflation stuff
%
% We have the following spaces.
%
% For the symmetries, we write:
%
% - P the symmetry group S(3) permuting the subsystems ABC, with P_AC the subgroup permuting only AC
%   We generate P using the three transpositions (A,B), (A,C), (B,C)
%
% - H the shift group S(2) shifting the ring by +3, or permuting the two triangles
%
% - L the symmetry group S(2) permuting the levels
%
% - D is the direct product P x H x L
%
% - D_AC is the direct product P_AC x H x L
%
% - T the torus group acting on the phases of rank 4 and dimension 6
%
% Triangle: C^8 three qubits A B C
% Symmetries: P+L+T, H has no action
%
% Two triangles: C^64 six qubits A1 B1 C1 A2 B2 C2
% Symmetries: P+H+L+T, where H permutes the two triangles
%
% Two triangles, PPT constraint: C^64
% Symmetries: P+L+T, H has no action
%
% Ring: C^64 six qubits A1 B1 C2 A2 B2 C2
% Symmetries: P+H+L+T, H acts by shifting the subsystems around the ring by three
%
% ACAC space: four qubits A1 C1 A2 C2
% Symmetries: P_AC+H+L+T *** we use a different group here!
%
% ABCB space: four qubits A1 B1 C1 B2 marginal of the ring
% Symmetries: P_AC+L+T, H has no action

% We define the groups

T6 = replab.T(6);
T6 = T6.withNames({'a0' 'b0' 'c0' 'a1' 'b1' 'c1'});
T = T6.subgroupWith('a0*b0*c0 = 1', 'a1*b1*c1 = 1');
P = replab.S(3);
P_AC = P.subgroup({[3 2 1]});
L = replab.S(2);
H = replab.S(2);

D = P.directProduct(H, L);

% Permutation of AB
gAB = {[2 1 3] [1 2] [1 2]};
actAB = T.automorphism('b0', 'a0', 'c0', 'b1', 'a1', 'c1');
% Permutation of AC
gAC = {[3 2 1] [1 2] [1 2]};
actAC = T.automorphism('c0', 'b0', 'a0', 'c1', 'b1', 'a1');
% Permutation of BC
gBC = {[1 3 2] [1 2] [1 2]};
actBC = T.automorphism('a0', 'c0', 'b0', 'a1', 'c1', 'b1');
% Shift
gH = {[1 2 3] [2 1] [1 2]};
actH = T.automorphism('a0', 'b0', 'c0', 'a1', 'b1', 'c1'); % identity
% Permutation of the two levels
gL = {[1 2 3] [1 2] [2 1]};
actL = T.automorphism('a1', 'b1', 'c1', 'a0', 'b0', 'c0');

D_AC = P_AC.directProduct(H, L);
D_AC_inj = D_AC.morphismByFunction(D, @(x) x); % we keep the same element

G = T.semidirectProductByFiniteGroup(D, 'preimages', {gAB, gAC, gBC, gH, gL}, 'images', {actAB, actAC, actBC, actH, actL});
G_AC = T.semidirectProductByFiniteGroup(D_AC, 'preimages', {gAC, gH, gL}, 'images', {actAC, actH, actL});
G_AC_inj = G_AC.morphismByFunction(G, @(x) x); % we keep the same element

% Triangle equivariant space
Trep_triangle = T.diagonalRepWith('a0 b0 c0', ...
                                  'a0 b0 c1', ...
                                  'a0 b1 c0', ...
                                  'a0 b1 c1', ...
                                  'a1 b0 c0', ...
                                  'a1 b0 c1', ...
                                  'a1 b1 c0', ...
                                  'a1 b1 c1');

imgAB_triangle = replab.Permutation.toMatrix([1 2 5 6 3 4 7 8]);
imgAC_triangle = replab.Permutation.toMatrix([1 5 3 7 2 6 4 8]);
imgBC_triangle = replab.Permutation.toMatrix([1 3 2 4 5 7 6 8]);
imgL_triangle = replab.Permutation.toMatrix([8 7 6 5 4 3 2 1]);
imgH_triangle = replab.Permutation.toMatrix([1 2 3 4 5 6 7 8]);

Drep_triangle = D.repByImages('C', 8, 'preimages', {gAB, gAC, gBC, gH, gL}, 'images', {imgAB_triangle, imgAC_triangle, imgBC_triangle, imgH_triangle, imgL_triangle});
rep_triangle = G.semidirectProductRep(Drep_triangle, Trep_triangle);
E_triangle = rep_triangle.hermitianInvariant;

% Two triangles equivariant space
Trep_twotriangles = kron(Trep_triangle, Trep_triangle);
imgAB_twotriangles = kron(imgAB_triangle, imgAB_triangle);
imgAC_twotriangles = kron(imgAC_triangle, imgAC_triangle);
imgBC_twotriangles = kron(imgBC_triangle, imgBC_triangle);
imgH_twotriangles = reshape(permute(reshape(eye(64), [8 8 64]), [2 1 3]), [64 64]);
imgL_twotriangles = kron(imgL_triangle, imgL_triangle);
Drep_twotriangles = D.repByImages('C', 64, 'preimages', {gAB, gAC, gBC, gH, gL}, 'images', {imgAB_twotriangles, imgAC_twotriangles, imgBC_twotriangles, imgH_twotriangles, imgL_twotriangles});
rep_twotriangles = G.semidirectProductRep(Drep_twotriangles, Trep_twotriangles);
E_twotriangles = rep_twotriangles.hermitianInvariant;

op_mapsTriangle = replab.equiop.generic(E_twotriangles, E_triangle, @(X) PartialTrace(X, [4 5 6], [2 2 2 2 2 2]));

% PPT equivariant space
Trep_ppt = kron(Trep_triangle, dual(Trep_triangle)); % the second copy has a partial transpose
Drep_ppt = D.repByImages('C', 64, 'preimages', {gAB, gAC, gBC, gH, gL}, 'images', {imgAB_twotriangles, imgAC_twotriangles, imgBC_twotriangles, eye(64), imgL_twotriangles});
rep_ppt = G.semidirectProductRep(Drep_ppt, Trep_ppt);
E_ppt = rep_ppt.hermitianInvariant;

op_ppt = replab.equiop.generic(E_twotriangles, E_ppt, @(X) PartialTranspose(X, [4 5 6], [2 2 2 2 2 2]), 'supportsSparse', true);

% ring equivariant space
Trep_ring = kron(Trep_triangle, Trep_triangle);
img_ring = @(p) reshape(permute(reshape(eye(64), [2 2 2 2 2 2 64]), [fliplr(7 - p) 7]), [64 64]);
imgAB_ring = img_ring([2 1 6 5 4 3]);
imgAC_ring = img_ring([6 5 4 3 2 1]);
imgBC_ring = img_ring([4 3 2 1 6 5]);
imgH_ring = img_ring([4 5 6 1 2 3]);
imgL_ring = kron(imgL_triangle, imgL_triangle);
Drep_ring = D.repByImages('C', 64, 'preimages', {gAB, gAC, gBC, gH, gL}, 'images', {imgAB_ring, imgAC_ring, imgBC_ring, imgH_ring, imgL_ring});
rep_ring = G.semidirectProductRep(Drep_ring, Trep_ring);
E_ring = rep_ring.hermitianInvariant;

% ACAC space, corresponds to A1C1A2C2 in the ring, A1C2A2C1
Trep_ACAC = T.diagonalRepWith('a0 c0 a0 c0', ...
                              'a0 c0 a0 c1', ...
                              'a0 c0 a1 c0', ...
                              'a0 c0 a1 c1', ...
                              'a0 c1 a0 c0', ...
                              'a0 c1 a0 c1', ...
                              'a0 c1 a1 c0', ...
                              'a0 c1 a1 c1', ...
                              'a1 c0 a0 c0', ...
                              'a1 c0 a0 c1', ...
                              'a1 c0 a1 c0', ...
                              'a1 c0 a1 c1', ...
                              'a1 c1 a0 c0', ...
                              'a1 c1 a0 c1', ...
                              'a1 c1 a1 c0', ...
                              'a1 c1 a1 c1');
img_ACAC = @(p) reshape(permute(reshape(eye(16), [2 2 2 2 16]), [fliplr(5 - p) 5]), [16 16]);
imgAC_ACAC = img_ACAC([4 3 2 1]);
imgH_ACAC = img_ACAC([3 4 1 2]);
imgL_ACAC = kron(kron([0 1; 1 0], [0 1; 1 0]), kron([0 1; 1 0], [0 1; 1 0]));
Drep_ACAC = D_AC.repByImages('C', 16, 'preimages', {gAC, gH, gL}, 'images', {imgAC_ACAC, imgH_ACAC, imgL_ACAC});
rep_ACAC = G_AC.semidirectProductRep(Drep_ACAC, Trep_ACAC);
E_ACAC = rep_ACAC.hermitianInvariant;

op_ACAC = replab.equiop.generic(E_ring, E_ACAC, @(X) PartialTrace(X, [2 5], [2 2 2 2 2 2]), 'sourceInjection', G_AC_inj);

% ABCB space, corresponds to A1B1C1B2 in the ring
Trep_ABCB = T.diagonalRepWith('a0 b0 c0 b0', ...
                              'a0 b0 c0 b1', ...
                              'a0 b0 c1 b0', ...
                              'a0 b0 c1 b1', ...
                              'a0 b1 c0 b0', ...
                              'a0 b1 c0 b1', ...
                              'a0 b1 c1 b0', ...
                              'a0 b1 c1 b1', ...
                              'a1 b0 c0 b0', ...
                              'a1 b0 c0 b1', ...
                              'a1 b0 c1 b0', ...
                              'a1 b0 c1 b1', ...
                              'a1 b1 c0 b0', ...
                              'a1 b1 c0 b1', ...
                              'a1 b1 c1 b0', ...
                              'a1 b1 c1 b1');
img_ABCB = @(p) reshape(permute(reshape(eye(16), [2 2 2 2 16]), [fliplr(5 - p) 5]), [16 16]);
imgAC_ABCB = img_ABCB([3 2 1 4]);
imgH_ABCB = eye(16);
imgL_ABCB = kron(kron([0 1; 1 0], [0 1; 1 0]), kron([0 1; 1 0], [0 1; 1 0]));
Drep_ABCB = D_AC.repByImages('C', 16, 'preimages', {gAC, gH, gL}, 'images', {imgAC_ABCB, imgH_ABCB, imgL_ABCB});
rep_ABCB = G_AC.semidirectProductRep(Drep_ABCB, Trep_ABCB);
E_ABCB = rep_ABCB.hermitianInvariant;

op_ABCB = replab.equiop.generic(E_ring, E_ABCB, @(X) PartialTrace(X, [4 6], [2 2 2 2 2 2]), 'sourceInjection', G_AC_inj);

% ABCB^T space, corresponds to A1B1C1B2 in the ring, with the B2 system partial transposed
Trep_ABCBT = T.diagonalRepWith('a0 b0 c0 b0^-1', ...
                               'a0 b0 c0 b1^-1', ...
                               'a0 b0 c1 b0^-1', ...
                               'a0 b0 c1 b1^-1', ...
                               'a0 b1 c0 b0^-1', ...
                               'a0 b1 c0 b1^-1', ...
                               'a0 b1 c1 b0^-1', ...
                               'a0 b1 c1 b1^-1', ...
                               'a1 b0 c0 b0^-1', ...
                               'a1 b0 c0 b1^-1', ...
                               'a1 b0 c1 b0^-1', ...
                               'a1 b0 c1 b1^-1', ...
                               'a1 b1 c0 b0^-1', ...
                               'a1 b1 c0 b1^-1', ...
                               'a1 b1 c1 b0^-1', ...
                               'a1 b1 c1 b1^-1');
rep_ABCBT = G_AC.semidirectProductRep(Drep_ABCB, Trep_ABCBT);
E_ABCBT = rep_ABCBT.hermitianInvariant;

op_ABCBT = replab.equiop.generic(E_ABCB, E_ABCBT, @(X) PartialTranspose(X, [4], [2 2 2 2]));