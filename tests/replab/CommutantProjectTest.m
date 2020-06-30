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
    W = S3.wreathProduct(replab.SignedPermutations.quaternionGroup);
    rep1 = W.primitiveRepFun(@(x) x.naturalRep);
    rep2 = rep1.decomposition;
    rep3 = rep2.asSimilarRep;
    X = randn(rep1.dimension, rep1.dimension);
    X1 = rep2.E_internal*rep1.commutant.project(rep2.B_internal*X*rep2.E_internal)*rep2.B_internal;
    X2 = rep2.commutant.project(X);
    X3 = rep3.commutant.project(X);
    assert(norm(X1 - X2) < replab.Parameters.doubleEigTol);
    assert(norm(X1 - X3) < replab.Parameters.doubleEigTol);
    assert(norm(X2 - X3) < replab.Parameters.doubleEigTol);
end

function test_complex_representations
    d = 10;
    S = replab.S(d);
    C = S.cyclicSubgroup;
    rep = C.naturalRep;
    rep1 = rep.decomposition;
    rep2 = rep1.asSimilarRep;
    X = randn(d, d);
    X1 = rep1.commutant.project(X);
    X2 = rep2.commutant.project(X);
    assert(norm(X1 - X2) < replab.Parameters.doubleEigTol);
end
