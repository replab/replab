function test_suite = IndexRelabelingRepTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    S3 = replab.S(3);
    rep = S3.indexRelabelingRep(2);
    test_suite = rep.laws.addTestCases(test_suite);
end
