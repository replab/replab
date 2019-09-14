function test_suite = nbVarsTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_cases
    global matrix231 matrix23451
    matrix1 = matrix23451;
    assert(length(matrix1.getVariables) == matrix1.nbVars);

    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    matrix2 = replab.CommutantVar.fromSdpMatrix(sdpvar(5,5,'hankel'), {[2 3 4 5 1]});
    assert(length(matrix2.getVariables) == matrix2.nbVars);
    assert(length(matrix2.getVariables) == matrix1.nbVars+5);
end
