function test_suite = RandomBagTest()
    disp(['Setting up tests in ', mfilename()]);
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
    Sn = replab.S(n);
    R = replab.bsgs.RandomBagWithImages(10, [g1' g2'], [], [], Sn, {g1 g2});
    for i = 1:100
        [x y] = R.sample;
        assert(isequal(x, y));
    end
end
