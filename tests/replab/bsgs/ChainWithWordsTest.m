function test_suite = ChainWithWordsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_leftCosetWord
    gens2 = {[2 3 4 5 1 6 7 8], [1 2 6 4 7 5 8 3]};
    G2 = replab.PermutationGroup.of(gens2{:});
    chain2 = replab.bsgs.ChainWithWords(G2, gens2);
    chain2.sgsWordQuick;
    chain2.setCompleted;

    v2 = [1 1 1 1 2 2 2 2];
    w2 = [1 2 1 1 2 1 2 2];
    P2 = G2.vectorFindPermutationsTo(v2, w2);
    a2long = chain2.word(P2.representative);
    [a2, b2] = chain2.wordCoset(P2);
    assert(length(a2long) > length(a2))
    assert(isequal(v2, w2(b2)))
end
