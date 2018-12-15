function test_suite = PermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_matrix()
    n = 100;
    assertTrue(isequal(replab.Perm.matrix(replab.Perm.identity(n)), eye(n)));
    for i = 1:100
        x = replab.Perm.random(n);
        y = replab.Perm.random(n);
        z = replab.Perm.compose(x, y);
        assertEqual(replab.Perm.matrix(z), replab.Perm.matrix(x) * replab.Perm.matrix(y));
    end
end

function test_group_laws()
    n = 100;
    cat = replab.cat.PermAsGroup(n);
    for i = 1:100
        cat.verifyLaws(@() randperm(n));
    end
end

function test_domain_action_laws()
    n = 100;
    cat = replab.cat.PermActingOnDomain(n);
    for i = 1:100
        cat.verifyLaws(@() randperm(n), @() randi(n));
    end
end

function test_pow()
    n = 100;
    for i = 1:100
        x = replab.Perm.random(n);
        m = randi(100);
        y = replab.Perm.identity(n);
        for j = 1:100
            y = replab.Perm.compose(y, x);
        end
        z = replab.Perm.pow(x, n);
        assertEqual(y, z);
    end
end

function test_isIdentity_identity()
    n = 100;
    assert(replab.Perm.isIdentity(replab.Perm.identity(100)));
    assert(~replab.Perm.isIdentity([1 3 2 4 5]));
end

function test_image()
    n = 100;
    for i = 1:100
        x = replab.Perm.random(n);
        y = replab.Perm.random(n);
        z = replab.Perm.compose(x, y);
        p = randi(n);
        imgy = replab.Perm.image(y, p);
        imgxy = replab.Perm.image(x, imgy);
        imgz = replab.Perm.image(z, p);
        assertEqual(imgxy, imgz);
    end
end

function test_vectorAction()
% test compatibility of action with group
    n = 100;
    for i = 1:100
        v = rand(n, 1);
        x = replab.Perm.random(n);
        y = replab.Perm.random(n);
        z = replab.Perm.compose(x, y);
        imgy = replab.Perm.vectorAction(y, v);
        imgxy = replab.Perm.vectorAction(x, imgy);
        imgz = replab.Perm.vectorAction(z, v);
        assertEqual(imgz, imgxy);
        Mz = replab.Perm.matrix(z);
        assertEqual(Mz * v, imgz);
    end
end

function test_selfAdjointMatrixAction()
    n = 10;
    for i = 1:100
        M = rand(n, n);
        x = replab.Perm.random(n);
        y = replab.Perm.random(n);
        z = replab.Perm.compose(x, y);
        imgy = replab.Perm.selfAdjointMatrixAction(y, M);
        imgxy = replab.Perm.selfAdjointMatrixAction(x, imgy);
        imgz = replab.Perm.selfAdjointMatrixAction(z, M);
        assertEqual(imgz, imgxy);
        rhoz = replab.Perm.matrix(z);
        assertEqual(rhoz * M * rhoz', imgz);
    end
end
