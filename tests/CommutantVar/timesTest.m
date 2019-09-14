function test_suite = timesTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    global matrix231 matrix23451
    matrix = matrix23451;
    difference = 2.*matrix - matrix - matrix;
    vars = [0 difference.getVariables];
    for j = 1:length(vars)
        coeffs = getBaseMatrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end

function test_inputs
    global matrix231 matrix23451
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix .* matrix);
end
