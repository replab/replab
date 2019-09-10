function test_suite = RandomBagTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_withPermutations
    n = 10;
    g1 = [2:n 1];
    g2 = [2 1 3:n];
    Sn = replab.Permutations(n);
    R = replab.bsgs1.RandomBag(10, [g1' g2'], [], [], Sn, {g1 g2});
    for i = 1:100
        [x y] = R.sample;
        assert(isequal(x, y));
    end
end
