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

% # Werner state symmetric extensions
%
% We symmetrize the computation in
%
% P. D. Johnson, "Compatible quantum correlations: Extension problems for Werner and isotropic states", PRA, vol. 88, no. 3, 2013
%
% https://journals.aps.org/pra/abstract/10.1103/PhysRevA.88.032323
%
%
%
%
% Initialize the RepLAB toolbox (be in the /replab directory or use `run path/replab/replab_init.m`)

% + tags=[]
replab_init
replab.globals.useReconstruction(1); % use new algorithms for decomposition
% -

%
% We declare the symmetry group $G$ as a direct product of those two groups:
%
% - $\mathcal{U}(2)$ acting on the subsystems
% - $\mathcal{S}(2)$ acting on the copies

U2 = replab.U(2);
S2 = replab.S(2);
G = U2.directProduct(S2)

% We now describe the action of the two factors on
% - the original state being tested
% - the symmetric extension we test the existence of

% +
% The action of U2 on the original two qubit state
U2_rep2 = kron(U2.definingRep, U2.definingRep);

% The action of U2 on the symmetric extension
U2_rep3 = kron(U2.definingRep, U2.definingRep, U2.definingRep);

% The action of S2 is trivial when there is a single copy of Bob
S2_rep2 = S2.trivialRep('C', 4);

% The action of S2 permutes the two copies of Bob (what we call a index relabeling)
S2_rep3 = kron(S2.trivialRep('C', 2), S2.indexRelabelingRep(2).complexification);
% -

% We construct the representations of the direct product using the representations above of the factor, which commute (important!).

rep2 = G.commutingFactorRepsRep('C', 4, {U2_rep2 S2_rep2});
rep3 = G.commutingFactorRepsRep('C', 8, {U2_rep3 S2_rep3});
H2 = rep2.hermitianInvariant;
H3 = rep3.hermitianInvariant;

% We define the partial trace operation (trust us).

ptFun = @(X) reshape(reshape(permute(reshape(X, [2 4 2 4]), [2 4 1 3]), [16 4])*[1; 0; 0; 1], [4 4]);
pt = replab.equiop.generic(H3, H2, ptFun);

% We define
%
% - the singlet state and noise as equivars (note that they are not variables)
% - the threshold `t` as a sdpvar
% - the symmetric extension we try to find `symExt`

singlet = replab.equivar(H2, 'value', [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]/2);
noise = replab.equivar(H2, 'value', eye(4)/4);
t = sdpvar;
rho = singlet*t + noise*(1-t);
symExt = replab.equivar(H3);

C = [pt(symExt) == rho
     issdp(symExt)
     issdp(rho)]
optimize(C, -t, sdpsettings('solver', 'sdpt3')) % force SDPT3 the default Octave solver has problems

double(t)
