function test_suite = blockTest()
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

function test_inputs
    generators = {[2 3 4 5 1]};
    matrix = replab.CommutantVar.fromPermutations(generators);
    shouldProduceAnError(@(x) matrix.block(-1));
    shouldProduceAnError(@(x) matrix.block(1.5));
    shouldProduceAnError(@(x) matrix.block(4));
end
