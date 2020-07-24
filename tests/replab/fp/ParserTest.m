function test_suite = ParserTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_parse_words
    assertEqual(replab.fp.parseLetters('[a,x]', {'a' 'x'}), [-1 -2 1 2]);
    assertEqual(replab.fp.parseLetters('a (x a)^2', {'a' 'x'}), [1 2 1 2 1]);
    assert(isempty(replab.fp.parseLetters('a a^-1', {'a' 'x'})));
end

function test_parse_correct_presentations
    [ok, names, relators] = replab.fp.Parser.parsePresentation('< a, x | a^2 = x^2 = 1,  a*x*a^-1 = x^-1 >');
    assert(ok);
    [ok, names, relators] = replab.fp.Parser.parsePresentation('< a, x >');
    assert(ok);
    [ok, names, relators] = replab.fp.Parser.parsePresentation('< a x >');
    assert(ok);
end
