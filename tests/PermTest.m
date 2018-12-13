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

function test_isIdentity()
    import replab.*
    n = 100;
    assert(Perm.isIdentity(Perm.identity(100)));
    assert(~Perm.isIdentity([1 3 2 4 5]));
end
