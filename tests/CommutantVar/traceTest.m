function test_suite = traceTest()
    try
        yalmip('version');
        try
            test_functions = localfunctions();
        catch
        end
        initTestSuite;
    catch
        warning('Yalmip not found in the path, some tests will be skipped');
        test_suite=MOxUnitTestSuite();
    end
end

function test_oneGroup
    % We do a sanity check with one group
    generators = {[2 3 4 5 1]};
    matrix = replab.CommutantVar.fromPermutations(generators);
    fullMatrix = matrix.fullMatrix;
    difference = trace(matrix) - trace(fullMatrix);
    vars = getvariables(difference);
    for j = 1:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end
