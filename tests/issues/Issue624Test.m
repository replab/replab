function test_suite = Issue624Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    G = replab.S(3);
    rep = G.naturalRep;
    space = rep.commutant;
    test_suite = space.laws.addTestCases(test_suite);
end
