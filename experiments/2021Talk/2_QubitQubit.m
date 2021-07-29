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

% # Qubit-qubit systems
%
% Initialize the RepLAB toolbox (be in the /replab directory or use `run path/replab/replab_init.m`)

replab_init
replab.globals.useReconstruction(1); % use new algorithms for decomposition

% Werner states and isotropic states are transformed by unitary matrices $\mathcal{U}(d)$.

d = 2;
Ud = replab.U(d)

% We take a random sample. A unitary group element is ... a unitary matrix. For groups in RepLAB, the random samples are actually taken from a uniform measure (the Haar measure).

Ud.sample

% If we represent group elements by ... unitary matrices, it's a quite straightforward representation called the "defining representation".

rep = Ud.definingRep

% This representation is irreducible, i.e. it doesn't have proper invariant subspaces.

rep.isIrreducible

% ## Werner states
%
% Werner states are invariant when applying the same unitary to the first and second subsystem.

rho_werner = kron(rep, rep)

% What are the invariant subspaces?

dec_werner = rho_werner.decomposition

dec_werner.basis

% ## Isotropic states
%
% Isotropic states are invariant when applying a unitary matrix to the first subsystem, and its conjugate to the second.

rho_iso = kron(rep, conj(rep))
dec_iso = rho_iso.decomposition
dec_iso.basis

% # Generalizations (exercices)
%
% ## Higher dimensions
%
% Both cases split in invariant subspaces of dimensions 1+3.
%
% What about dimension 3? dimension 4?
%
% ## Cloning machines
%
% The Choi state of an isotropic cloning machine would transform as:
%
% `kron(rep, rep, ..., rep, conj(rep), conj(rep), ..., conj(rep))`
%
% where the number of `rep` is the number of output copies, and the number of `conj(rep)` is the number of input copies.
%
% In 2018, the decomposition was made analytically for the case of a single `conj(rep)` copy: http://stacks.iop.org/1751-8121/51/i=12/a=125202
