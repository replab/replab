function test_suite = IrreducibleTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_irreducible_commutant
    S20 = replab.S(20);
    rep = S20.naturalRep;
    dec = rep.decomposition;
    assertEqual(dec.irrep(1,1).dimension, 1);
    C = dec.commutant;
    X = randn(20, 20);
    X1 = C.projectAndReduceFromParent(X);
end
