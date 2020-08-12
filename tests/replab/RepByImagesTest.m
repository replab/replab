function test_suite = RepByImagesTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    % cyclotomics
    G = replab.AbstractGroup.parsePresentation('<a, x | a^8 = x^2 = 1, x a x^-1 = a^-1 >');
    sqrt2 = replab.cyclotomic.sqrtRational(2);
    img_a = [1/sqrt2 -1/sqrt2
             1/sqrt2 1/sqrt2];
    img_x = [1 0; 0 -1];
    rep = G.repByImages('R', 2, {img_a img_x});
    test_suite = rep.laws.addTestCases(test_suite);
end

function test_symbolic
    S8 = replab.S(8); % we write D_16 as a subgroup of S_8
    a = [2 3 4 5 6 7 8 1]; % cyclic rotation
    x = [8 7 6 5 4 3 2 1]; % reflection
    D16 = S8.subgroup({a x});
    s2 = sqrt(sym(2));
    sigma_a = [1 -1; 1 1]/s2;
    sigma_x = [1 0; 0 -1];
    sigma = D16.repByImages('R', 2, {sigma_a sigma_x});
    img1 = sigma.image(S8.compose(a, x)); % should be [1 1; 1 -1]/sqrt(2)
    img2 = [1 1; 1 -1]/s2;
    assert(norm(double(img1) - double(img2)) == 0);
end

function test_stays_sparse
    S3 = replab.PermutationGroup.of([2 3 1], [2 1 3]);
    M1 = sparse([0 0 1; 1 0 0; 0 1 0]);
    M2 = sparse([0 1 0; 1 0 0; 0 0 1]);
    rep = S3.repByImages('R', 3, {M1 M2});
    for i = 1:5
        I = rep.image_internal(S3.sample);
        assert(issparse(I));
    end
end