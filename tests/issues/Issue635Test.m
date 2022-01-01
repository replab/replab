function test_suite = Issue635Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_issue()
    S3 = replab.S(3);
    rep = S3.signRep;
    G = rep.kernel; % G should be C(3)
    assertFalse(G.isTrivial); % fails
end
