function test_suite = valueTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    % We do some sanity checks
    global matrix231 matrix23451
    matrix = matrix23451;
    
    assert(norm(isnan(value(matrix))-1) == 0);
end
