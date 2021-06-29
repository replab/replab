function test_suite = StabilizersTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_subgroup
    if ReplabTestParameters.onlyFastTests
        n = 1;
        d = 3;
    else
        n = 10;
        d = 10;
    end
    for i = 1:n
        S = replab.S(d);
        G = S.randomProperSubgroup(2);
        g = G.sample;
        P = replab.Partition.fromVector(randi(3, 1, d));
        O = G.orderedPartitionStabilizer(P);
        U = G.unorderedPartitionStabilizer(P);
        assert(O.isSubgroupOf(U));
    end
end
