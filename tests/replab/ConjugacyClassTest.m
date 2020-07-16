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
    classes = replab.nfg.conjugacyClassesByOrbits(G);
    for i = 1:length(classes)
        cc1 = classes{i};
        cc2 = replab.ConjugacyClass(G, cc1(:,1)');
        assert(size(cc1, 2) == cc2.size);
    end
end
