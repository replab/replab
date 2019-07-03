Q = replab.SignedPermutations.quaternionGroup;
S3 = replab.Permutations(3);
W = replab.WreathProductGroup(S3, Q);
rep = W.primitiveRep(Q.naturalRep);
sub = replab.rep1.decompose(rep);
S = sub{1}    
W = replab.rep1.enforceQuaternionEncoding(S);
W' * S.sample * W
