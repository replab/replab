function test_suite = nbVarsTest()
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
    assert(length(matrix.getVariables) == matrix.nbVars);
end
