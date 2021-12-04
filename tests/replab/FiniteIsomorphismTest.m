function test_suite = FiniteIsomorphismTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S2 = replab.PermutationGroup.of([2 1]);
    SS2 = replab.SignedPermutationGroup.signedSymmetric(1);
    iso = S2.isomorphismByImages(SS2, 'images', {-1});
    test_Suite = iso.laws.addTestCases(test_suite);
end
