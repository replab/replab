% Declare the symmetry group: U(2) acting on each subsystem
% We do not use permutation of subsystems as it complicates the PPT constraint

G = replab.U(2);
% The representation of the two qubit state
rep = kron(G.definingRep, G.definingRep);
% The representation acting on the partially transposed state
repT = kron(G.definingRep, dual(G.definingRep));

% Spaces of equivariant Hermitian matrices
H = rep.hermitianInvariant;
HT = repT.hermitianInvariant;

% The partial transpose linear map
pptFun = @(X) reshape(permute(reshape(full(X), [2 2 2 2]), [3 2 1 4]), [4 4]);
% Upgraded as an equivariant super operator from the space H to the space HT
ppt = replab.equiop(H, HT, pptFun);

% Define the singlet state as an equivar
singlet = replab.equivar(H, 'value', [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]/2);
% Define the noise as an equivar
noise = replab.equivar(H, 'value', eye(4)/4);

% Standard sdpvar
t = sdpvar;
% equivar*sdpvar works, but not sdpvar*equivar (as sdpvar would provide the * op and doesn't handle equivar)
rho = singlet*t + noise*(1-t);

% we use the syntax sdp(EV) to define a semidefinite positive constraint, where EV is an equivar
% we use the syntax sdpvar(EV) to recover the sdpvar in the original (i.e. non symmetry adapted) basis

C = [real(trace(sdpvar(rho))) == 1
     sdp(rho)
     sdp(ppt(rho))];

% note that this is a linear program now
optimize(C, -t)
