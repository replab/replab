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

% # Werner state separability
%
% Initialize the RepLAB toolbox (be in the /replab directory or use `run path/replab/replab_init.m`)

replab_init
replab.globals.useReconstruction(1); % use new algorithms for decomposition

% Declare the symmetry group: $G = \mathcal{U}(2)$ acting on each subsystem. We do not implement the permutation of subsystems as it complicates the PPT constraint.

G = replab.U(2);

% This is the representation leaving the Werner state invariant.

rep = kron(G.definingRep, G.definingRep);

% Now, the partially transposed state has a "transpose" on the second subsystem; i.e. we take the conjugate representation.

repT = kron(G.definingRep, conj(G.definingRep));

% We define spaces of Hermitian matrices invariant under the representations. Those are spaces of matrices that transform as $X \rightarrow \rho(g) X \rho(g)^\dagger$ for each of the representations we defined.

H = rep.hermitianInvariant
HT = repT.hermitianInvariant

% We define the partial transpose linear map. We could have used http://www.qetlab.com/PartialTranspose as well. Note that we permute the indices 1,3 as Matlab's reshape does not follow the kron convention (long story).

ptFun = @(X) reshape(permute(reshape(X, [2 2 2 2]), [3 2 1 4]), [4 4]);

% We tell RepLAB that the partial transpose is a super operator from the space H to the space HT; by doing so, RepLAB knows that the operator is compatible with the action of the group (`H` and `HT` must have been defined using representations of the same group to use this syntax).

pt = replab.equiop.generic(H, HT, ptFun)

% Define the singlet state as an invariant matrix (equivariant variable in the math jargon).

singlet = replab.equivar(H, 'value', [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]/2)

% Same for the noise

noise = replab.equivar(H, 'value', eye(4)/4)

% The visibility is a standard sdpvar from YALMIP.

t = sdpvar;

% Below, the syntax `equivar*sdpvar` works, but not `sdpvar*equivar`, which is why the scalar is on the right (as `sdpvar` would provide the `*` operator and doesn't handle equivar)

rho = singlet*t + noise*(1-t);

% We use the syntax `sdp(X)` to define a semidefinite positive constraint, where `X` is an equivar.

C = [issdp(rho)
     issdp(pt(rho))]

% note that this is a linear program now

optimize(C, -t, sdpsettings('solver', 'sdpt3')) % force SDPT3 the default Octave solver has problems

% The separability threshold is:

double(t)


