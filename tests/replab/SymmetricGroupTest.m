function test_suite = SymmetricGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_conjugacy_classes
    G = replab.S(10);
    C1 = G.conjugacyClasses;
    H = replab.PermutationGroup.of([2:10 1], [2 1 3:10]);
    C2 = H.conjugacyClasses;
    for i = 1:C1.nClasses
        G.assertEqv(C1.classes{i}.representative, C2.classes{i}.representative);
        assert(C1.classes{i}.representativeCentralizer == C2.classes{i}.representativeCentralizer);
    end
end
