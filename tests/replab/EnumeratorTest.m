function test_suite = EnumeratorTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_display
    G = replab.Permutations(60);
    description = replab.shortStr(G.elements);
    assert(~isempty(description));
end
