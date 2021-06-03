% ---
% jupyter:
%   jupytext:
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.11.2
%   kernelspec:
%     display_name: Octave
%     language: octave
%     name: octave
% ---

% # Introduction: groups and representations
%
% Initialize the RepLAB toolbox (be in the /replab directory or use `run path/replab/replab_init.m`)

replab_init

% Construct the symmetric group.

S3 = replab.S(3)

% Let's take a random element from the group; it's a permutation row vector

S3.sample

% Construct the representation that acts on $\vec{P} = (P(1), P(2), P(3))$ in two different ways: the first one is a shortcut, the second one uses explicit images. We check the soundness of our construction.

% +
rho = S3.naturalRep;
rho = S3.repByImages('R', 3, 'preimages', {[2 3 1] [2 1 3]}, 'images', {[0 0 1; 1 0 0; 0 1 0] [0 1 0; 1 0 0; 0 0 1]})
rho.check % Run automated tests

% this one is wrong
%rho = S3.repByImages('R', 3, 'preimages', {[2 3 1] [2 1 3]}, 'images', {[0 1 0; 0 0 1; 1 0 0] [0 1 0; 1 0 0; 0 0 1]})
%rho.errorBound
% -

% Let's take the image of a group element. It's the corresponding permutation matrix.

g = [3 2 1];
rho.image(g)

% Finally, let us look at the invariant subspaces.

S2 = replab.S(2)
rho = S2.naturalRep
dec = rho.decomposition

% The representation has an invariant subspace `[1,1,1]` and the subspace `[2 -1 -1; 0 1 -1]` is invariant as well.

% ## Manipulating representations

% Now, imagine the same group is acting on $P(a_1, a_2) \in \mathbb{R}^9$ which represents two successive outcomes of the box. We could construct `rho2` by computing the explicit images, but the tensor product of the representation works as well here.

rho2 = kron(rho, rho, rho) % tensor product

dec2 = rho2.decomposition

% Let's look at the discovered basis. The first two columns correspond to the two copies of the trivial representation. Those are invariant vectors. The third column is also an invariant subspace. The remaining columns decompose in a more complex way (three copies of an irreducible representation of dimension 2, the standard representation of $\mathcal{S}_3$). 

dec2.basis


