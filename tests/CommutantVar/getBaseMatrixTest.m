function test_suite = getBaseMatrixTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231 + eye(3);
    vars = [0 matrix.getVariables];
    assert(length(vars) == 3);
    assert(max(max(abs(matrix.getBaseMatrix(vars(1)) - eye(3)))) < replab.globals.doubleEigTol);
    assert(max(max(abs(matrix.getBaseMatrix(vars(2)) - 1/3))) < replab.globals.doubleEigTol);
    assert(max(max(abs(matrix.getBaseMatrix(vars(3)) + 1/3 - eye(3)))) < replab.globals.doubleEigTol);
end
