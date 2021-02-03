S3 = replab.S(3);
W = S3.wreathProduct(replab.QuaternionGroup());
rep = W.primitiveRepFun(@(G,i) G.naturalRep);
subs = rep.complexSplit;
