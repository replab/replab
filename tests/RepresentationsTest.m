function test_suite = RepresentationsTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.naturalRepresentation;
    I = rho.irreducible;
    I1 = rho.irreducible.component(1);
    I2 = rho.irreducible.component(2);
    test_suite = replab.RealRepLaws(rho).addTestCases(test_suite);
    test_suite = replab.RealRepLaws(I).addTestCases(test_suite);
    test_suite = replab.RealRepLaws(I1).addTestCases(test_suite);
    test_suite = replab.RealRepLaws(I2).addTestCases(test_suite);
    G = replab.Permutations(5);
    test_suite = replab.RepLaws(G.standardRep).addTestCases(test_suite);
end

function test_symmetric_group_representations
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.naturalRepresentation;
    I = rho.irreducible;
    assertEqual(I.nComponents, 2);
    assertEqual(cellfun(@(c) c.dimension, I.components), [1 3]);
    assertEqual(cellfun(@(c) c.multiplicity, I.components), [1 1]);
end
