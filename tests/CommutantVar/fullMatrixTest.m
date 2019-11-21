function test_suite = fullMatrixTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    M = sdpvar(2);
    matrix = replab.CommutantVar.fromSdpMatrix(M, {[2 1]});
    assert(norm(matrix.fullMatrix-M) == 0);
end
