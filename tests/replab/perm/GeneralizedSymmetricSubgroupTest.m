function test_suite = GeneralizedSymmetricSubgroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    G = replab.perm.GeneralizedSymmetricGroupType(3, 3);
    G = G.parentGroup;
    test_suite = G.laws.addTestCases(test_suite);
    H = G.randomProperSubgroup(2);
    test_suite = H.laws.addTestCases(test_suite);
end
