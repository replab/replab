function test_suite = SymmetricGroupIrrepsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ~replab.init.cyclolab().works
        return
    end

    if ReplabTestParameters.onlyFastTests
        youngs = {[2 1]};
    else
        youngs = {[1] [2] [1 1] [3] [2 1] [1 1 1] [4]};% [3 1] [2 2] [2 1 1] [1 1 1 1]}; % last part commented out as per issue #247
    end

    for i = 1:length(youngs)
        Y = youngs{i};
        Sn = replab.S(sum(Y));
        % This part of the test is currently commented for most i as of issue #247
        rep = Sn.irrep(Y, 'specht');
        test_suite = rep.laws.addTestCases(test_suite);
        rep = Sn.irrep(Y, 'seminormal');
        test_suite = rep.laws.addTestCases(test_suite);
        rep = Sn.irrep(Y, 'orthogonal');
        test_suite = rep.laws.addTestCases(test_suite);
    end
end
