function test_suite = Issue451Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_bug
% Before the fix this would take 15! property tests
    n = 15;
    S = replab.S(n);
    G = S.vectorStabilizer(1:n);
    assert(G.order == 1);
end
