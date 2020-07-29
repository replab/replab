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
    G = Sn.randomProperSubgroup(2);
    classes = replab.perm.conjugacyClassesByOrbits(G);
    for i = 1:length(classes)
        cc1 = classes{i};
        cc2 = replab.ConjugacyClass(G, cc1(:,1)');
        assert(size(cc1, 2) == cc2.size);
    end
end

function test_conjugacyClass_number
    d = [5 6 7 8 9];
    n = [7 11 15 22 30]; % partition function
    if ReplabTestParameters.onlyFastTests
        d = d(1:3);
        n = n(1:3);
    end
    for i = 1:length(d)
        S = replab.S(d(i));
        C = S.conjugacyClasses;
        sz = sum(cellfun(@(c) double(c.size), C));
        assert(sz == S.order);
        ni = length(S.conjugacyClasses);
        assert(n(i) == ni);
    end
end
