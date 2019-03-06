function test_suite = DivisionAlgebraTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_display
    G = replab.Permutations(60);
    description = G.elements.str;
    assert(~isempty(description));
end
