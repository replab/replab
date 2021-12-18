function test_suite = FiniteMorphismTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.PermutationGroup.of([2 3 1], [2 1 3]);
    S2 = replab.S(2);
    signMorphism1 = S3.morphismByImages(S2, 'images', {[1 2] [2 1]});
    test_suite = signMorphism1.laws.addTestCases(test_suite);

    SS2 = replab.SignedPermutationGroup.signedSymmetric(1);
    signMorphism2 = S3.morphismByImages(SS2, 'images', {1 -1});
    test_suite = signMorphism2.laws.addTestCases(test_suite);
end
