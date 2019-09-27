function test_suite = SedumiDataTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_CHSH
    thisPath = mfilename('fullpath');
    [testsPath, ~, ~] = fileparts(thisPath);
    dataPath = fullfile(testsPath, 'CHSH_sedumi.mat');
    S = replab.SedumiData.fromMatFile(dataPath);
    [F h x] = S.toYalmip;
    solvesdp(F, h);
    obj = -2*sqrt(2);
    assert(abs(double(h) - obj) < 1e-6);
end
