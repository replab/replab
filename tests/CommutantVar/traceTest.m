function test_suite = traceTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_oneGroup
    % We do a sanity check with one group
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
    fullMatrix = matrix.fullMatrix;
    difference = trace(matrix) - trace(fullMatrix);
    vars = getvariables(difference);
    for j = 0:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end
