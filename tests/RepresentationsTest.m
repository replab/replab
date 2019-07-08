function test_suite = RepresentationsTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    % test orbit decomposition
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.naturalRep;
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
    % test standard representation
    G = replab.Permutations(5);
    test_suite = replab.RepLaws(G.standardRep).addTestCases(test_suite);
    % test quaternion type representations
    Q = replab.SignedPermutations.quaternionGroup;
    S3 = replab.S(3);
    W = replab.WreathProductGroup(S3, Q);
    rho = W.primitiveRep(Q.naturalRep);
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
    rho = W.imprimitiveRep(Q.naturalRep);
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
end

function test_symmetric_group_representations
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.naturalRep;
    I = rho.decomposition;
    assertEqual(I.nComponents, 2);
    assertEqual(cellfun(@(c) c.copyDimension, I.components), [1 3]);
    assertEqual(cellfun(@(c) c.multiplicity, I.components), [1 1]);
end

function test_representation_of_cyclic_group
    C12 = replab.Permutations(12).cyclicSubgroup;
    rep = C12.naturalRep;
    d = cellfun(@(iso) iso.copyDimension, rep.decomposition.components);
    m = cellfun(@(iso) iso.multiplicity, rep.decomposition.components);
    t = cellfun(@(iso) iso.copy(1).realDivisionAlgebra.d, rep.decomposition.components);
    assertEqual(d, [1 1 2 2 2 2 2]);
    assertEqual(m, [1 1 1 1 1 1 1]);
    assertEqual(t, [1 1 2 2 2 2 2]);
end
