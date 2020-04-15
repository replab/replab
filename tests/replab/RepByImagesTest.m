function test_suite = RepByImagesTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_symbolic
    S8 = replab.S(8); % we write D_16 as a subgroup of S_8
    a = [2 3 4 5 6 7 8 1]; % cyclic rotation
    x = [8 7 6 5 4 3 2 1]; % reflection
    D16 = S8.subgroup({a x});
    s2 = sqrt(sym(2));
    sigma_a = [1 -1; 1 1]/s2;
    sigma_x = [1 0; 0 -1];
    sigma = D16.repByImages('R', 2, {sigma_a sigma_x}, {sigma_a' sigma_x'});
    img1 = sigma.image(S8.compose(a, x)); % should be [1 1; 1 -1]/sqrt(2)
    img2 = [1 1; 1 -1]/s2;
    assert(double(simplify(norm(img1-img2))) == 0);
end
