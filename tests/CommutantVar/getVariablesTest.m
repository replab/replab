function test_suite = getVariablesTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    global matrix231 matrix23451
    matrix = matrix23451;
    list1 = matrix.getVariables;
    list2 = getvariables(matrix.fullMatrix);
    assert(isequal(sort(list1), sort(list2)));
end
