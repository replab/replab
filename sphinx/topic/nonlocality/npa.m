% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:light
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.3.2
%   kernelspec:
%     display_name: Octave
%     language: octave
%     name: octave
% ---

% # Symmetric NPA relaxations
%
% We initialize the RepLAB library.

addpath([pwd, '/../../..']);
replab_init('verbose', 0);

% We examine an upper bound on the quantum maximum of the CHSH inequality. We use the NPA hierarchy, and assume the moment matrix has been constructed in the
% symmetric subspace, as in [arXiv:1808.09598](https://arxiv.org/abs/1808.09598).
%
% We first declare the one and only variable, before constructing the moment matrix.

y = sdpvar; % y is actually y_A0B0, but we save space by identifying y with it
C = [ 1  0  0  0  0
      0  1  0  0  0
      0  0  1  0  0
      0  0  0  1  0
      0  0  0  0  1];
A = [ 0  0  0  0  0
      0  0  0  1  1
      0  0  0  1 -1
      0  1  1  0  0
      0  1 -1  0  0];
X = C + A*y;
I_CHSH = 4*y; % it is <A0B0> + <A0B1> + <A1B0> - <A1B1>, and y_A1B1 = -y_A0B0.
optimize(X >= 0, -I_CHSH, sdpsettings('verbose', 0)); % sign change to maximize, and we don't show solver output
double(I_CHSH)

% To help with symmetrization, we'll use an equivalent form of the constraint X, and check the result.

% Now, we want to find a change of basis that brings X into a block-diagonal form. This is easier to check on the basis matrices `C` and `A`, i.e. we want to find `U` such that `U*C*U'` and `U*A*U'` are block-diagonal. Let us use RepLAB for that.
% First, remark that the moment matrix is invariant under the following signed permutations:
%
% $\vec{v} = \left( 1, A_0, A_1, B_0, B_1 \right) \rightarrow  \left( 1, B_0, B_1, A_0, A_1 \right)$
%
% $\vec{v} = \left( 1, A_0, A_1, B_0, B_1 \right) \rightarrow  \left( 1, -A_0, -A_1, -B_0, -B_1 \right)$
%
% $\vec{v} = \left( 1, A_0, A_1, B_0, B_1 \right) \rightarrow  \left( 1, A_1, A_0, B_0, -B_1 \right)$
%
% We now show how to use RepLAB to find the change of basis for the symmetry group that includes only the first symmetry.

g1 = [1 4 5 2 3]; % permutation of parties
g2 = [1 -2 -3 -4 -5]; % sign flip everywhere, note the signed permutation convention
g3 = [1 3 2 4 -5]; % additional symmetry

% We build the symmetry group from those generators. For the permutation of parties only, the order (=size) of the group is 2. When adding the other symmetries, the order of the group should be 16.

nElements = 5;G = replab.SignedPermutationGroup.of(g1, g2, g3);
G.order

% If necessary, we can also expand the group elements (use `G.elements.at(2)` to get the second element for example).

G.elements

% Now, the representation we need is composed of signed permutation matrices; as we are considering a group of signed permutations, this is the group *natural representation*. Representations in RepLAB are described using the images of the generators. By calling `G.sample`, we get random elements from the group, and by calling `rep.image(g)` we get the image of a group element.

rep = G.naturalRep

g = G.sample

% Now, we decompose the representation into irreducible components. For that, we call `rep.decomposition`. In the answer I(m)xR(d), we express that the component has $m$ copies (multiplicity) of a `R`eal representation of dimension `d`. We then play with the indices of the component, and copy.
%
% Hint: explore `rep.decomposition.component(i).irrep(j)` for $(i,j) = (1,1), (2,1), (3,1)$.

rep.decomposition

rep.decomposition.component(1).irrep(1)

% We now ask for the change of basis matrix, and verify that `A` and `C` block-diagonalize.

U = rep.decomposition.basis
U'*C*U
U'*A*U

% Once we have found the change of basis matrix, we transform our moment matrix `X`, and naturally it should have a block diagonal structure following the one of `C` and `A` -- for that, it may be necessary to kill small off-block coefficients.

X = (U'*C*U) + (U'*A*U)*y; % equivalent to the constraint above
optimize(X >= 0, -I_CHSH); % sign change to maximize
double(I_CHSH)

% As an exercice, recover an algebraic basis from `U`, guess the form of the semidefinite program when fully block-diagonal (i.e. for the full group of CHSH symmetries).
% **Remark that now the problem can be solved by hand to recover the** $2\sqrt{2}$ **bound, as the SDP matrix is fully diagonal, and thus corresponds to linear inequalities.**
