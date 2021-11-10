function test_suite = realeigTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
% function test_commutant_sample
%     G = replab.PermutationGroup.cyclic(7);
%     rep = G.naturalRep;
%     tol = 1e-6;
%     for i = 1:10
%         X = rep.commutant.sample;
%         [V,D,W] = replab.numerical.realeig(X);
%         assert(norm(X*V - V*D, 'fro') <= tol);
%         assert(norm(W'*X - D*W', 'fro') <= tol);
%         D = diag(D);
%         assert(all(D(2:end) >= D(1:end-1)));
%     end
% end
function test_normal_samples
    tol = 1e-6;
    for i = 1:10
        X = randn(9,9);
        [V,D,W] = replab.numerical.realeig(X);
        assert(norm(X*V - V*D, 'fro') <= tol);
        assert(norm(W'*X - D*W', 'fro') <= tol);
        D = diag(D);
        assert(all(D(2:end) >= D(1:end-1)));
    end
end
