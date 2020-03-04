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
    X1 = I.asSimilarRep.commutant.project(X);
    X2 = I.commutant.project(X1);
    assert(norm(X1 - X2) < replab.Parameters.doubleEigTol);
end

function test_complex_representations
    d = 10;
    S = replab.S(d);
    C = S.cyclicSubgroup;
    rep = C.naturalRep;
    rep1 = rep.decomposition;
    rep2 = rep1.asSimilarRep;
    X = randn(d, d);
    X1 = rep1.commutant.project(X)
    X2 = rep2.commutant.project(X);
    assert(norm(X1 - X2) < replab.Parameters.doubleEigTol);
end
