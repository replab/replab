function test_suite = Issue414Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_bug
    S5 = replab.S(5);
    G1 = S5.subgroup({[2 3 1 4 5], [2 3 4 5 1]});
    G2 = S5.subgroup({[2 3 1 4 5], [1 2 4 5 3]});
    assert(G1.recognize.atlasEntry.group.order == 60);
    assert(G2.recognize.atlasEntry.group.order == 60);
end
