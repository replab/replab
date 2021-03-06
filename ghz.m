n = 3;
d = 3;
Sn = replab.S(n);
Sd = replab.S(d);
% full torus is reshape of phases(n, d)
% torus is reshape of phases(n-1, d)
torusToFull = [eye(n-1); -ones(1, n-1)];
fullToTorus = [eye(n-1) zeros(n-1, 1)];
Sn_images = cell(1, Sn.nGenerators);
Sn_nr = Sn.naturalRep;
for i = 1:Sn.nGenerators
    Sn_images{i} = fullToTorus*Sn_nr.image(Sn.generator(i))*torusToFull;
end
Sn_rep = Sn.repByImages('R', n-1, 'preimages', Sn.generators, 'images', Sn_images);
G = Sd.directProduct(Sn);
torusRep = G.tensorFactorRep('R', {Sd.naturalRep Sn_rep});
ghzGroup = replab.TorusGroup.semidirectProductFromRep(torusRep);
stateFiniteRep = G.commutingProductFactorRep('R', d^n, {Sd.naturalRep.tensorPower(n), Sn.indexRelabelingRep(d)});
tm = zeros(d^n, n-1, d);
for i = 1:n-1
    tm = reshape(tm, [d^(i-1) d d^(n-i) n-1 d]);
    for j = 1:d
        tm(:,j,:,i,j) = tm(:,j,:,i,j) + 1;
    end
end
tm = reshape(tm, [d^(n-1) d n-1 d]);
for j = 1:d
    tm(:,j,:,j) = tm(:,j,:,j) - 1;
end
tm = reshape(tm, [d^n (n-1)*d]);
contRep = ghzGroup.N.mapRep(tm);
rep = ghzGroup.productRep(stateFiniteRep.complexification, contRep);