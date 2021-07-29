function test_suite = IrreducibleTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_irreducible_commutant
    if ReplabTestParameters.onlyFastTests
        d = 3;
    else
        d = 20;
    end
    Sd = replab.S(d);
    rep = Sd.naturalRep;
    dec = rep.decomposition;
    assertEqual(dec.irrep(1,1).dimension, 1);
    C = dec.commutant;
    X = randn(d, d);
    X1 = C.project(X);
    % TODO: add test
end

function test_import_export
    G = replab.S(3);
    parent = G.naturalRep;
    dec = parent.decomposition('exact').squeeze;
    data = dec.export;
    dec1 = replab.Irreducible.import(parent, data);
    dec1.check;
    dec = parent.decomposition;
    data = dec.export;
    dec1 = replab.Irreducible.import(parent, data);
    dec1.check;
end