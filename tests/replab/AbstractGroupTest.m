function test_suite = AbstractGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        allNs = [10];
    else
        allNs = [0 1 2 5 10];
    end

    G = replab.AbstractGroup.fromPresentation('< x | x^3 = 1 >');
    test_suite = G.laws.addTestCases(test_suite);
    G = replab.AbstractGroup.fromPresentation('< x, a | a^3 = x^2 = a x a x = 1 >');
    test_suite = G.laws.addTestCases(test_suite);
end

function test_identity_generator
    assertExceptionThrown(@() replab.AbstractGroup.make({'x'}, {'x'}));
end