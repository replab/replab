function test_suite = RepByImagesTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    % cyclotomics
    G = replab.AbstractGroup.parsePresentation('<a, x | a^8 = x^2 = 1, x a x^-1 = a^-1 >');
    sqrt1_2 = replab.cyclotomic.sqrtRational(1,2);
    img_a = [sqrt1_2 -sqrt1_2
             sqrt1_2 sqrt1_2];
    img_x = replab.cyclotomic.fromDoubles([1 0; 0 -1]);
    rep = G.repByImages('R', 2, {img_a img_x});
    test_suite = rep.laws.addTestCases(test_suite);
end

function test_cyclotomic
    S8 = replab.S(8); % we write D_16 as a subgroup of S_8
    a = [2 3 4 5 6 7 8 1]; % cyclic rotation
    x = [8 7 6 5 4 3 2 1]; % reflection
    D16 = S8.subgroup({a x});
    s2 = replab.cyclotomic.sqrtRational(1, 2);
    sigma_a = replab.cyclotomic.sqrtRational(1, 2) * [1 -1; 1 1];
    sigma_x = replab.cyclotomic.fromDoubles([1 0; 0 -1]);
    sigma = D16.repByImages('R', 2, 'preimages', {a x}, 'images', {sigma_a sigma_x});
    img1 = sigma.image(S8.compose(a, x)); % should be [1 1; 1 -1]/sqrt(2)
    img2 = replab.cyclotomic.sqrtRational(1, 2) * [1 1; 1 -1];
    assert(norm(double(img1) - double(img2)) == 0);
end

function test_stays_sparse
    S3 = replab.PermutationGroup.of([2 3 1], [2 1 3]);
    M1 = sparse([0 0 1; 1 0 0; 0 1 0]);
    M2 = sparse([0 1 0; 1 0 0; 0 0 1]);
    rep = S3.repByImages('R', 3, {M1 M2});
    for i = 1:5
        I = rep.image(S3.sample, 'double/sparse');
        assert(issparse(I));
    end
end

function test_pauli_group_is_unitary_1_design
% Pauli group from presentation
    G = replab.AbstractGroup.parsePresentation('< x,y,z,i | i^4 = x^2 = y^2 = z^2 = 1, x i = i x, y i = i y, z i = i z, x y = i z, y z = i x, z x = i y>');
    images = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1], 1i*[1 0; 0 1]};
    rep = G.repByImages('C', 2, 'preimages', {'x' 'y' 'z' 'i'}, 'images', images);
    assert(rep.decomposition.nComponents == 1);
end

function test_pauli_group_is_not_unitary_2_design
% Pauli group from presentation
    G = replab.AbstractGroup.parsePresentation('< x,y,z,i | i^4 = x^2 = y^2 = z^2 = 1, x i = i x, y i = i y, z i = i z, x y = i z, y z = i x, z x = i y>');
    images = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1], 1i*[1 0; 0 1]};
    images2 = cellfun(@(I) kron(I, I), images, 'uniform', 0);
    rep = G.repByImages('C', 4, 'preimages', {'x' 'y' 'z' 'i'}, 'images', images2);
    assert(rep.decomposition.nComponents == 4);
end
