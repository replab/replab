function test_suite = SymmetricGroupIrrepsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        youngs = {[2 1]};
    else
        youngs = {[1] [2] [1 1] [3] [2 1] [1 1 1] [4]};% [3 1] [2 2] [2 1 1] [1 1 1 1]}; % last part commented out as per issue #247
    end

    if exist('syms')
        for i = 1:length(youngs)
            Y = youngs{i};
            if i == 1
                % This part of the test is currently commented for most i as of issue #247
                rep = replab.sym.specht(Y);
                test_suite = rep.laws.addTestCases(test_suite);
                rep = replab.sym.youngSeminormal(Y);
                test_suite = rep.laws.addTestCases(test_suite);
            end
            rep = replab.sym.youngOrthogonal(Y);
            test_suite = rep.laws.addTestCases(test_suite);
        end
    end
end
