function test_suite = CosetsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    D8 = replab.PermutationGroup.dihedral(8);
    C8 = replab.PermutationGroup.cyclic(8);
    leftCosets = D8/C8;
    rightCosets = C8\D8;
    test_suite = leftCosets.laws.addTestCases(test_suite);
    test_suite = rightCosets.laws.addTestCases(test_suite);

    SS3 = replab.SignedPermutationGroup.signedSymmetric(3);
    sub = SS3.subgroup({[2 3 1], [-1 -2 -3]});
    leftCosets = SS3/sub;
    rightCosets = sub\SS3;
    test_suite = leftCosets.laws.addTestCases(test_suite);
    test_suite = rightCosets.laws.addTestCases(test_suite);

    S5 = replab.S(5);
    H = S5.randomProperSubgroup;
    K = S5.randomProperSubgroup;
    g = S5.sample;
    doubleCoset = H * g * K;
    test_suite = doubleCoset.laws.addTestCases(test_suite);
end
