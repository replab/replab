function test_suite = lmtimesTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    difference = 2*matrix - matrix - matrix;
    vars = [0 difference.getVariables];
    for j = 1:length(vars)
        coeffs = getBaseMatrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Parameters.doubleEigTol);
    end
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix * matrix);
end
