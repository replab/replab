function test_suite = PermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_matrix()
    import replab.*
    n = 100;
    assertTrue(isequal(Perm.matrix(Perm.identity(n)), eye(n)));
    for i = 1:100
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.compose(x, y);
        assertEqual(Perm.matrix(z), Perm.matrix(x) * Perm.matrix(y));
    end
end

function test_compose()
    import replab.*
    n = 100;
    for i = 1:100
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.random(n);
        xy = Perm.compose(x, y);
        yz = Perm.compose(y, z);
        assertEqual(Perm.compose(x, yz), Perm.compose(xy, z));
    end
end

function test_isIdentity_identity()
    import replab.*
    n = 100;
    assert(Perm.isIdentity(Perm.identity(100)));
    assert(~Perm.isIdentity([1 3 2 4 5]));
end

function test_inverse()
    import replab.*
    n = 100;
    for i = 1:100
        x = Perm.random(n);
        y = Perm.random(n);
        assertEqual(Perm.identity(n), Perm.compose(x, Perm.inverse(x)));
        xy = Perm.compose(x, y);
        yIxI = Perm.compose(Perm.inverse(y), Perm.inverse(x));
        assertEqual(Perm.inverse(xy), yIxI);
    end
end

function test_image()
    import replab.*
    n = 100;
    for i = 1:100
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.compose(x, y);
        p = randi(n);
        imgy = Perm.image(y, p);
        imgxy = Perm.image(x, imgy);
        imgz = Perm.image(z, p);
        assertEqual(imgxy, imgz);
    end
end

function test_vectorAction()
% test compatibility of action with group
    import replab.*
    n = 100;
    for i = 1:100
        v = rand(n, 1);
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.compose(x, y);
        imgy = Perm.vectorAction(y, v);
        imgxy = Perm.vectorAction(x, imgy);
        imgz = Perm.vectorAction(z, v);
        assertEqual(imgz, imgxy);
        Mz = Perm.matrix(z);
        assertEqual(Mz * v, imgz);
    end
end

function test_matrixAction()
    import replab.*
    n = 10;
    for i = 1:100
        M = rand(n, n);
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.compose(x, y);
        imgy = Perm.matrixAction(y, M);
        imgxy = Perm.matrixAction(x, imgy);
        imgz = Perm.matrixAction(z, M);
        assertEqual(imgz, imgxy);
        rhoz = Perm.matrix(z);
        assertEqual(rhoz * M * rhoz', imgz);
    end
end