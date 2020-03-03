function test_suite = MorphismsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.Permutations(3);
    S27 = replab.Permutations(27);
    phi = S3.indexRelabelingMorphism(3);
    laws = replab.GroupMorphismLaws(phi, S3, S27);
    test_suite = laws.addTestCases(test_suite);
end
