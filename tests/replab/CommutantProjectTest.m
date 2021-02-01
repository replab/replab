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
    W = S3.wreathProduct(replab.QuaternionGroup);
    rep1 = W.primitiveRepFun(@(x) x.naturalRep);
    rep2 = rep1.decomposition;
    rep3 = rep2.asSimilarRep;
    X = randn(rep1.dimension, rep1.dimension);
    X1 = rep2.projection*rep1.commutant.project(rep2.injection*X*rep2.projection)*rep2.injection;
    X2 = rep2.commutant.project(X);
    X3 = rep3.commutant.project(X);
    assert(norm(X1 - X2) < replab.globals.doubleEigTol);
    assert(norm(X1 - X3) < replab.globals.doubleEigTol);
    assert(norm(X2 - X3) < replab.globals.doubleEigTol);
end

function test_complex_representations
    d = 10;
    S = replab.S(d);
    C = replab.CyclicGroup(d);
    rep = C.naturalRep;
    rep1 = rep.decomposition;
    rep2 = rep1.asSimilarRep;
    X = randn(d, d);
    X1 = rep1.commutant.project(X);
    X2 = rep2.commutant.project(X);
    assert(norm(X1 - X2) < replab.globals.doubleEigTol);
end
