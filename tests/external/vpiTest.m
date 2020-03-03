function test_suite = vpiTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_vpiLongNumbers
    if ReplabTestParameters.onlyFastTests
        return;
    end

    replab.longStr(replab.Permutations(100).elements);
    % used to throw because num2str of @vpi split numbers into lines
end
