function sub = decompose(rep)
% Decomposes the given representation into irreducible subrepresentations
    assert(isa(rep, 'replab.Rep'));
    sub1 = replab.rep1.orbitDecomposition(rep);
    trivial = cell(1, 0);
    nontrivial = cell(1, 0);
    for i = 1:length(sub1)
        s1 = sub1{i};
        [Ut Ur] = replab.rep1.extractTrivial(s1);
        for j = 1:size(Ut, 2)
            trivial{1, end+1} = s1.subRep(Ut(:, j)).in(rep);
        end
        s2 = s1.subRep(Ur);
        sub3 = replab.rep1.decomposeUsingCommutant(s2);
        for j = 1:length(sub3)
            s3 = sub3{j};
            nontrivial{1, end+1} = s3.in(rep);
        end
    end
    sub = horzcat(trivial, nontrivial);
end
