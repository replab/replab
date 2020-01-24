function test_suite = gtTest()
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
    
    % We do some sanity checks
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]}, 'symmetric', 'real');
    warn = evalc('assert(length(matrix > 2) == 3)');
    if replab.compat.isOctave
        assert(isequal(warn(1:9), 'warning: '))
    else
        assert(isequal(warn(3:11), 'Warning: '))
    end
end
