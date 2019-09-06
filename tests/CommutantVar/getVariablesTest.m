function test_suite = getVariablesTest()
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
    list1 = matrix.getVariables;
    list2 = getvariables(matrix.fullMatrix);
    assert(isequal(sort(list1), sort(list2)));
end
