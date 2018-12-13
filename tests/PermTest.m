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
    for i = 1:100
        x = Perm.random(n);
        y = Perm.random(n);
        z = Perm.compose(x, y);
        assertEqual(Perm.matrix(z), Perm.matrix(x) * Perm.matrix(y));
    end
end
