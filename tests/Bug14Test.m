function test_suite = Bug14Test()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    group = replab.SignedPermutations(4).subgroup({[1 -3 2 4]});
    I = group.naturalRep.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
end
