C = replab.PermutationGroup.cyclic(3);
rep = C.naturalRep;
ctC = C.characterTable('C');
ctR = C.characterTable('R');
els = C.elements.toCell;
ch1 = ctC.character(2);
ch2 = ctC.character(3);
ir = ctR.irreps{2};
P1 = zeros(2, 2);
P2 = zeros(2, 2);
for i = 1:length(els);
    g = els{i};
    P1 = P1 + conj(ch1.value(g))*ir.image(g,'exact');
    P2 = P2 + conj(ch2.value(g))*ir.image(g,'exact');
end
P1 = P1 / length(els);
P2 = P2 / length(els);
P1 = zeros(2, 2);
P2 = zeros(2, 2);
for i = 1:length(els);
    g = els{i};
    P1 = P1 + conj(ch1.value(g))*ir.image(g,'exact');
    P2 = P2 + conj(ch2.value(g))*ir.image(g,'exact');
end
P1 = P1 / length(els);
P2 = P2 / length(els);
T1 = zeros(6, 6);
T2 = zeros(6, 6);
for i = 1:length(els)
    g = els{i};
    T1 = T1 + kron(P1*ir.image(g,'exact')*P1, rep.image(g,'exact'));
    T2 = T2 + kron(P2*ir.image(g,'exact')*P2, rep.image(g,'exact'));
end
T1 = T1 / length(els);
T2 = T2 / length(els);
T1 = permute(reshape(T1, [3 2 3 2]), [1 3 2 4]);
T2 = permute(reshape(T2, [3 2 3 2]), [1 3 2 4]);
