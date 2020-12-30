function test_suite = SubRepTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    F = [ 2  0
         -1  1
         -1 -1];
    G = [2 -1 -1
         0  3 -3]/6;
    S3 = replab.S(3);
    rep = S3.naturalRep.subRep(F, 'projection', G);
    test_suite = rep.laws.addTestCases(test_suite);
end
