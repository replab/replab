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
    X1 = C.projectAndReduceFromParent(X);
    % TODO: add test
end
