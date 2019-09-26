function test_suite = CommutantProjectTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_quaternion_representations
    S3 = replab.S(3);
    W = S3.wreathProduct(replab.signed.Permutations.quaternionGroup);
    rep = W.primitiveRepFun(@(x) x.naturalRep);
    X = randn(rep.dimension, rep.dimension);
    I = rep.decomposition;
    X1 = I.asConjugateRep.commutant.project(X); 
    X2 = I.asRep.commutant.project(X1);
    assert(norm(X1 - X2) < replab.Settings.doubleEigTol);
end

function test_complex_representations
    S20 = replab.S(20);
    C20 = S20.cyclicSubgroup;
    rep = C20.naturalRep;
    X = randn(rep.dimension, rep.dimension);
    I = rep.decomposition;
    X1 = I.asConjugateRep.commutant.project(X); 
    X2 = I.asRep.commutant.project(X1);
    assert(norm(X1 - X2) < replab.Settings.doubleEigTol);
end
