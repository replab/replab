% # Tutorial for physicists
%
% The unitary group $U(n)$ describes the possible change of basis in the
% complex Hilbert space of dimension $n$. For $n=2$, this captures the
% possible choices of basis in the qubit space $C^2$. For a system composed
% of two qubits, we can similarly define the effect of a joint change of
% basis performed on both subsystems simultaneously. It is known that a
% single state is invariant under such joint change of basis. Here, we
% identify this state by extracting the subspace of $(C^2)^{\otimes 2}$
% which is invariant under the joint change of basis for both subsystems.

% Before trying any of the *RepLAB* commands, we must initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## The unitary group representation
%
% Changes of bases for one system are described by the group $U(2)$

d = 2;
U2 = replab.U(d);

% We construct the defining representation of this group, which acts on $C^2$:

U2Rep = U2.definingRep;

% ## Tensor product of two representations
%
% We can now construct the representation which acts jointly on two
% subsystems of dimension 2:

U2TensorRep = kron(U2Rep, U2Rep);

% To identify the subspaces which are invariant under this group, we decompose the representation:

dec = U2TensorRep.decomposition.nice;

% The decomposition has

dec.nComponents

% components, of dimension

dec.component(1).irrepDimension

% and

dec.component(2).irrepDimension

% These are the antisymmetric and symmetric subspaces respectively. The change of basis into the first component
% identifies the antisymetric subspace, also known as the *singlet state*:

singletBasis = dec.component(1).basis
rest = dec.component(2).basis
