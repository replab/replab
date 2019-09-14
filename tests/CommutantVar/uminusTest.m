function test_suite = uminusTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    % We do some sanity checks
    global matrix231 matrix23451
    matrix = matrix23451;
    opposite = -matrix;
    difference = matrix + opposite;
    vars = [0 difference.getVariables];
    for j = 1:length(vars)
        coeffs = getBaseMatrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end
