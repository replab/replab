function test_suite = ParserTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_parse_correct_presentations
    P = replab.fp.Parser;
    [ok, names, relators] = P.parsePresentation('< a, x | a^2 = x^2 = 1,  a*x*a^-1 = x^-1 >');
    assert(ok);
    [ok, names, relators] = P.parsePresentation('< a, x >');
    assert(ok);
    [ok, names, relators] = P.parsePresentation('< a x >');
    assert(ok);
end
