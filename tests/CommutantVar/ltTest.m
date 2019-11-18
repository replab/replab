function test_suite = ltTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]}, 'symmetric', 'real');
    assert(length(matrix < 2) == 3);
end
