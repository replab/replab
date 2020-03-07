function test_suite = IsotypicTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_changeOfBasis
    S3 = replab.S(3);
    rep = S3.naturalRep.tensorPower(2);
    rep1 = rep.decomposition;
    c3 = rep1.component(3);
    t1 = rep.subRep([1 0 0 0 1 0 0 0 1]');
    t1.isIrreducible = true;
    t2 = rep.subRep([0 1 1 1 0 1 1 1 0]');
    t2.isIrreducible = true;
    s = rep.subRep([0 1 -1 -1 0 1 1 -1 0]');
    s.isIrreducible = true;
    d1 = rep.subRep([2 0 0 0 -1 0 0 0 -1; -1 0 0 0 2 0 0 0 -1]');
    d1.isIrreducible = true;
    d2 = rep.subRep([0 -1 1 -1 0 1 0 0 0; 0 1 -1 0 0 0 -1 1 0]');
    d2.isIrreducible = true;
    d3 = rep.subRep([0 2 2 -1 0 -1 -1 -1 0; 0 -1 -1 2 0 2 -1 -1 0]');
    d3.isIrreducible = true;
    iso1 = replab.Isotypic.fromIrreps(rep, {t1 t2});
    iso3 = replab.Isotypic.fromIrreps(rep, {d1 d2 d3});
    replab.IsotypicLaws(iso1).check;
    replab.SubRepLaws(s).check;
    replab.IsotypicLaws(iso3).check;
    Miso = replab.domain.Matrices(iso3.field, iso3.dimension, iso3.dimension);
    M = replab.domain.Matrices(d1.field, d1.dimension, d1.dimension);
    for i = 1:3
        for j = 1:3
            [A Ainv] = iso3.changeOfBasis(i, j);
            iso3h = iso3.changeIrrepBasis(j, A, Ainv);
            Miso.assertEqv(iso3h.E_internal * iso3h.B_internal, eye(iso3h.dimension));
            for k = 1:S3.nGenerators
                g = S3.generator(k);
                M.assertEqv(A*iso3.irrep(j).image(g)*Ainv, iso3.irrep(i).image(g));
                M.assertEqv(iso3h.irrep(j).image(g), iso3.irrep(i).image(g));
            end
        end
    end
end
