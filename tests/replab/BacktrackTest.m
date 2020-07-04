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
    for j = 1:3
        G = replab.nfg.randomSubgroup(G);
    end
    v = randi(3, 1, n);
    s = G.sample;
    w = [];
    w(s) = v;
    g = G.findPermutationTo(v, w);
    assert(~isempty(g));
    assert(isequal(v(g), w));
end
