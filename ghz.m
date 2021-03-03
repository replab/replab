n = 3;
d = 2;
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
