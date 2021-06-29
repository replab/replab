function test_suite = minusTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    difference = matrix - matrix;
    vars = [0 difference.getVariables];
    for j = 1:length(vars)
        coeffs = getBaseMatrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.globals.doubleEigTol);
    end
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    shouldProduceAnError(@(x) matrix - rand(5));
    shouldProduceAnError(@(x) matrix - rand(3));
    shouldProduceAnError(@(x) rand(5) - matrix);
    shouldProduceAnError(@(x) rand(3) - matrix);
end
