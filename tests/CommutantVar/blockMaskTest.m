function test_suite = blockMaskTest()
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
