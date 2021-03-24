
G = replab.U(2);
rep = kron(G.definingRep, G.definingRep);
repT = kron(G.definingRep, G.definingRep.');

H = rep.hermitianInvariant;
HT = repT.hermitianInvariant;

singlet = replab.equivar(H, 'value', [0 0 0 0; 0 1 -1 0; 0 -1 1 0; 0 0 0 0]/2);
noise = replab.equivar(H, 'value', eye(4)/4);

t = sdpvar;

rho = singlet*t + noise*(1-t);

pptFun = @(X) reshape(permute(reshape(full(X), [2 2 2 2]), [3 2 1 4]), [4 4]);
ppt = replab.equiop(H, HT, pptFun);

C = [real(trace(sdpvar(rho))) == 1
     sdp(rho)
     sdp(ppt(rho))];
optimize(C, -t)
