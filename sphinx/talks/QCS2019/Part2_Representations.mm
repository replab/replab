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

% # 2. Constructing representations with RepLAB
%
% This is part II of the companion notebook to the RepLAB talk at the [Quantum Causal Structures](http://www.cs.ox.ac.uk/conferences/QCS2019/) workshop.

run ../../../replab_init.m % Init RepLAB library

% ## Two measurement settings with two outcomes
% We come back to our example for conditional probability distributions for two measurements settings with two outcomes.

piI = [3 4 1 2];
piO1 = [2 1 3 4];
GAlice = replab.S(4).subgroup({piI piO1})

rhoI = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
rhoO1 = [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
rho = GAlice.repByImages('R', 4, {rhoI rhoO1}, {rhoI' rhoO1'})

% We now verify that `rho` is a legitimate representation.

replab.RepLaws(rho).check % verify that this defines a proper representation

% We also define and verify the other representation `sigma`.

sigmaI = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
sigmaO1 = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
sigma = GAlice.repByImages('R', 4, {sigmaI sigmaO1}, {sigmaI' sigmaO1'})

% ## Combining representations
% We illustrate the examples of states invariant under $U(2)$. At the time of writing this tutorial, RepLAB did support compact groups, only finite groups; as of March 2020, the support of continuous groups is still experimental. Thus, we use the Clifford group as an approximation of $U(2)$. The approximation stays valid for tensor products of up to three copies of the basic representation, as the Clifford group is a 3-design for $U(2)$.
%
% The Clifford group is described using a permutation representation obtained with [GAP System](https://www.gap-system.org) on 48 points.
%
% We first get the (abstract) Clifford group on a single qubit and its representation using 2x2 complex matrices:

[clifford cliffordRep] = replab.quantum.clifford_qudit(2)

% We now investigate tensor products of `cliffordRep`, as an approximation of $\rho \otimes \rho$ and later $\rho \otimes \overline{\rho}$.

% ### The singlet and two qubit Werner states
% We investigate qubit-qubit states invariant under $\rho \otimes \rho$, where $\rho$ is the standard representation of $U(2)$.

UxU = kron(cliffordRep, cliffordRep)

UxU.decomposition.nice

% Without surprise, we recover the singlet and triplet (dimensions 1 and 3). Let's investigate.

UxU.decomposition.nice.component(1).irrep(1)

UxU.decomposition.nice.component(2).irrep(1)

% ### Choi of depolarizing channels / isotropic states
% We investigate qubit-qubit states invariant under $\rho \otimes \overline{\rho}$, where $\rho$ is the standard representation of $U(2)$.

UxconjU = kron(cliffordRep, conj(cliffordRep))

UxconjU.decomposition.nice

UxconjU.decomposition.nice.component(1).irrep(1)

UxconjU.decomposition.nice.component(2).irrep(1)
