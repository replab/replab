function test_suite = timesTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451H;
    difference = 2.*matrix - matrix - matrix;
    vars = [0 difference.getVariables];
    for j = 1:length(vars)
        coeffs = getBaseMatrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Parameters.doubleEigTol);
    end
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix .* matrix);
    shouldProduceAnError(@(x) matrix .* ones(1,3));
end
