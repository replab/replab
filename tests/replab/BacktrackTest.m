function test_suite = BacktrackTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_backtrack
    n = 10;
    G = replab.S(n);
    G = G.randomProperSubgroup(3);
    s = randi(3, 1, n);
    g = G.sample;
    t(g) = s;
    g1 = G.vectorFindPermutationsTo(s, t).representative;
    assert(~isempty(g));
    assert(isequal(t(g), s));
end
