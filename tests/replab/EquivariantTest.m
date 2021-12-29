function test_suite = EquivariantTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    G = replab.S(3);
    rep1 = G.naturalRep;
    E = rep1.antilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.bilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.commutant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.hermitianInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.symmetricInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.sesquilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    rep1 = rep1.complexification;

    E = rep1.antilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.bilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.commutant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.hermitianInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.symmetricInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.sesquilinearInvariant;
    test_suite = E.laws.addTestCases(test_suite);
    if ~replab.init.cyclolab().works
        return
    end
    rep1 = G.naturalRep;
    rep2 = G.irrep([2 1], 'specht');
    E = rep1.equivariantTo(rep2);
    test_suite = E.laws.addTestCases(test_suite);
    E = rep1.equivariantFrom(rep2);
    test_suite = E.laws.addTestCases(test_suite);
end
