function test_suite = PermutationGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_pointwiseStabilizer_with_trivial_result
% Used to crash when the result was trivial
    n = 5;
    S = replab.S(n);
    G = S.pointwiseStabilizer(1:n);
    assert(G.order == 1);
end
