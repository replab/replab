%% The antisymmetric subspace of U(2)xU(2) and the singlet state
%
% The unitary group $U(n)$ describes the possible change of basis in the
% complex Hilbert space of dimension $n$. For $n=2$, this captures the
% possible choices of basis in the qubit space $C^2$. Using *RepLAB*'s
% decomposition algorithm, we identify here the subspace of the two-qubits
% Hilbert space $(C^2)^{\otimes 2}$ which is invariant under the join
% change of basis for both subsystems.

%%
% Before trying any of the *RepLAB* commands, we must initialize the library:
replab_init

%% The unitary group representation
% Changes of bases for one system are described by the group $U(2)$
d = 2;
U2 = replab.UnitaryGroup(d);
%%
% We construct the defining representation of this group, which acts on
% $C^2$:
U2Rep = U2.definingRep;

%% Tensor product of two representations
% We can now construct the representation which acts jointly on two
% subsystems of dimension 2:
U2TensorRep = kron(U2Rep, U2Rep);
%%
% To identify the subspaces which are invariant under this group, we
% decompose the representation:
dec = U2TensorRep.decomposition;
%%
% The decomposition has
dec.nComponents
%%
% components, of dimension
dec.component(1).irrepDimension
%%
% and
dec.component(2).irrepDimension
%%
% These are the antisymmetric and symmetric subspaces respectively. The
% change of basis into the first component identifies the antisymetric
% subspace, also known as the *Singlet state*:
phiMinus = dec.component(1).U
