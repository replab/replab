function test_suite = blockTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_inputs
    global matrix231 matrix23451
    matrix = matrix23451;
    shouldProduceAnError(@(x) matrix.block(-1));

    if ReplabTestParameters.onlyFastTests
        return;
    end
    
    shouldProduceAnError(@(x) matrix.block(1.5));
    shouldProduceAnError(@(x) matrix.block(4));
end
