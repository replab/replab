function test_suite = vpiTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_vpiLongNumbers
    replab.longStr(replab.Permutations(100).elements)
    % used to throw because num2str of @vpi split numbers into lines
end
