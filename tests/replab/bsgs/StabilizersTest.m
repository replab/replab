function test_suite = StabilizersTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_subgroup
    n = 10;
    d = 10;
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
