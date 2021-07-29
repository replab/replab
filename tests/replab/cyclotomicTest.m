function test_suite = cyclotomicTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_parse_print
    c = replab.cyclotomic({'E(4) + E(7)'});
    s = strtrim(num2str(c));
    c1 = replab.cyclotomic({s});
    assert(c == c1);
end

function test_parse_print_rand
    for i = 1:10
        a = replab.cyclotomic.rand;
        b = replab.cyclotomic(num2str(a));
        assert(a == b);
    end
end
