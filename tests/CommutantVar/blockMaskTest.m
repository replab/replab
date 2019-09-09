function test_suite = blockMaskTest()
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
    % We do some sanity checks
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
    mask = matrix.blockMask;
    assert(issparse(mask));
    assert(nnz(mask) == 9);
    
    matrix = replab.CommutantVar.fromPermutations({[1 3 2]});
    mask = matrix.blockMask;
    assert(issparse(mask));
    assert(nnz(mask) == 5);
end
