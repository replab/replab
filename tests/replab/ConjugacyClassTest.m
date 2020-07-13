function test_suite = ConjugacyClassTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_centralizer_order
    n = 8;
    Sn = replab.S(n);
    G = replab.nfg.randomSubgroup(Sn);
    G = replab.nfg.randomSubgroup(G);
    C = G.conjugacyClasses;
    for i = 1:length(C)
        g = C{i}.representative;
        ctr = G.centralizer(C{i}.representative);
        assert(C{i}.cardinality * ctr.order == G.order);
    end
end
