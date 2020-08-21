function test_suite = Issue400Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_bug
    S3 = replab.S(3);
    sub = S3.subgroup({});
    sub.recognize; % should not throw
end
