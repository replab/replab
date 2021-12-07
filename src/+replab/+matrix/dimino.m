function H = dimino(generators, identity)
    d = size(identity, 1);
    H = replab.numerical.CycloSet(d, d);
    H.insert(identity);
    for i = 1:length(generators)
        H = replab.matrix.inductiveStep(H, generators(1:i), identity);
    end
end
