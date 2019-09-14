function test_suite = traceTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_oneGroup
    % We do a sanity check with one group
    global matrix231 matrix23451
    matrix = matrix23451;
    fullMatrix = matrix.fullMatrix;
    difference = trace(matrix) - trace(fullMatrix);
    vars = [0 getvariables(difference)];
    for j = 1:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end
