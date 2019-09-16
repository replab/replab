function test_suite = blockMaskTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    mask = matrix.blockMask;
    assert(issparse(mask));
    assert(nnz(mask) == 9);
    
    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    matrix = matrix231;
    mask = matrix.blockMask;
    assert(issparse(mask));
    assert(nnz(mask) == 5);
end
