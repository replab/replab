function test_suite = rdivideTest()
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

function test_cases
    generators = {[2 3 4 5 1]};
    matrix = replab.CommutantVar.fromPermutations(generators);

    difference = (2*matrix)./2 - matrix;
    vars = difference.getVariables;
    for j = 1:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end

function test_inputs
    matrix = replab.CommutantVar.fromPermutations({[3 2 1]});
    shouldProduceAnError(@(x) matrix./matrix);
end
