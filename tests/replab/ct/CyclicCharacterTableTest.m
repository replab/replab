function test_suite = CyclicCharacterTableTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    if ~replab.init.cyclolab().works
        return
    end
    G = replab.atl.Cyclic.make(3);
    test_suite = G.realCharacterTable.laws.addTestCases(test_suite);
    test_suite = G.complexCharacterTable.laws.addTestCases(test_suite);
end
