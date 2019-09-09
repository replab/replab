function test_suite = sizeTest()
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

function test_oneGroup
    % We do a sanity check with one group
    matrix = replab.CommutantVar.fromPermutations({[2 3 4 5 1]});
    s12 = size(matrix);
    assert(s12(1) == size(matrix,1));
    assert(s12(2) == size(matrix,2));
end

function test_inputs
    matrix = replab.CommutantVar.fromPermutations({[3 2 1]});
    shouldProduceAnError(@(x) size(matrix, 3));
    shouldProduceAnError(@(x) size(1, matrix));
end
