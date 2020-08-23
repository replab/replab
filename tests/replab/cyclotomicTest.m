function test_suite = cyclotomicTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_parse_print
    c = replab.cyclotomic.fromStrings({'E(4) + E(7)'});
    s = strtrim(num2str(c));
    c1 = replab.cyclotomic.fromStrings({s});
    assert(c == c1);
end
