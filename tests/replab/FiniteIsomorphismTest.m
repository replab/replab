function test_suite = FiniteIsomorphismTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S2 = replab.PermutationGroup.of([2 1]);
    SS2 = replab.SignedSymmetricGroup(1);
    iso = S2.isomorphismByImages(SS2, {-1});
    test_Suite = iso.laws.addTestCases(test_suite);
end
