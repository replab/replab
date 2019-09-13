function test_suite = gtTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    % We do some sanity checks
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
    assert(length(matrix > 2) == 3);
end
