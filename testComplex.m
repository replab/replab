C3 = replab.Permutations(3).cyclicGroup;
S3 = replab.Permutations(3);
W = replab.WreathProductGroup(S3, C3);
rep = W.primitiveRep(C3.standardRep);
sub = replab.rep1.decompose(rep);
if sub{1}.dimension == 2
    S = sub{2};
else
    S = sub{1};
end
    
W = replab.rep1.enforceComplexEncoding(S);
W' * S.sample * W