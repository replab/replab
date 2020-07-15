function test_suite = SetTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_sort
    G = replab.DihedralGroup(10);
    s = replab.perm.Set(G.domainSize);
    s.insert(G.chain.allElements);
    s.sort;
    s.check
end
