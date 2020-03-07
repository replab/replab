function test_suite = Issue68Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_bug
    if ReplabTestParameters.onlyFastTests
        return;
    end

    group = replab.signed.Permutations(12).subgroup({[1 3 2 5 6 7 4 9 8 12 10 11]});
    I = group.naturalRep.decomposition;
    assert(isa(I, 'replab.Irreducible'));
end
