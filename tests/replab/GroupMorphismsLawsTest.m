function test_suite = GroupMorphismsLawsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.S(3);
    S27 = replab.S(27);
    phi = S3.indexRelabelingMorphism(3);
    laws = replab.laws.GroupMorphismLaws(phi, S3, S27);
    test_suite = laws.addTestCases(test_suite);
end
