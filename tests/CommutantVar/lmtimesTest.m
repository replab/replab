function test_suite = lmtimesTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});

    difference = 2*matrix - matrix - matrix;
    vars = difference.getVariables;
    for j = 0:length(vars)
        coeffs = getbasematrix(difference, vars(j));
        assert(norm(coeffs(:)) <= replab.Settings.doubleEigTol);
    end
end

function test_inputs
    matrix = replab.CommutantVar.fromPermutations({[3 2 1]});
    shouldProduceAnError(@(x) matrix * matrix);
end
