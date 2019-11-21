function test_suite = traceTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    
    matrix = matrix23451;
    fullMatrix = matrix.fullMatrix;
    difference = trace(matrix) - trace(fullMatrix);
    vars = [0 getvariables(difference)];
    for j = 1:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Parameters.doubleEigTol);
    end
    
    matrix = matrix23451H;
    fullMatrix = matrix.fullMatrix;
    difference = trace(matrix) - trace(fullMatrix);
    vars = [0 getvariables(difference)];
    for j = 1:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Parameters.doubleEigTol);
    end
    
end
