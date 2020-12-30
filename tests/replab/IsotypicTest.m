function test_suite = IsotypicTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.S(3);
    rep = S3.naturalRep.tensorPower(2);
    dec = rep.decomposition;
    c3 = dec.component(3);
    t1 = rep.subRep([1 0 0 0 1 0 0 0 1]', 'isIrreducible', true, 'trivialDimension', 1);
    t2 = rep.subRep([0 1 1 1 0 1 1 1 0]', 'isIrreducible', true, 'trivialDimension', 1);
    s = rep.subRep([0 1 -1 -1 0 1 1 -1 0]', 'isIrreducible', true, 'trivialDimension', 0);
    d1 = rep.subRep([2 0 0 0 -1 0 0 0 -1; -1 0 0 0 2 0 0 0 -1]', 'trivialDimension', 0, 'isIrreducible', true);
    d2 = rep.subRep([0 -1 1 -1 0 1 0 0 0; 0 1 -1 0 0 0 -1 1 0]', 'trivialDimension', 0, 'isIrreducible', true);
    d3 = rep.subRep([0 2 2 -1 0 -1 -1 -1 0; 0 -1 -1 2 0 2 -1 -1 0]', 'trivialDimension', 0, 'isIrreducible', true);
    iso1 = replab.Isotypic.fromIrreps(rep, {t1 t2}, 1, true);
    iso2 = replab.Isotypic.fromIrreps(rep, {s}, 1, true);
    iso3 = replab.Isotypic.fromIrreps(rep, {d1 d2 d3}, 2, false);
    test_suite = iso1.laws.addTestCases(test_suite);
    test_suite = iso2.laws.addTestCases(test_suite);
    test_suite = iso3.laws.addTestCases(test_suite);
    test_suite = iso3.harmonize.laws.addTestCases(test_suite);
end
