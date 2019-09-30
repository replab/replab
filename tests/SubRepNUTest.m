function test_suite = SubRepNUTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    F = [2 -1 -1
         0  1 -1];
    G = [ 2  0
         -1  3
         -1 -3]/6;
    S3 = replab.S(3);
    rep = S3.definingRep.subRep(F, G);
    test_suite = replab.SubRepNULaws(rep).addTestCases(test_suite);
end
