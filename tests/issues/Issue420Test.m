function test_suite = Issue420Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_bug
    S3 = replab.S(3);
    G = S3.directProduct(S3, S3);
    assert(~isequal(G, {}));
end
