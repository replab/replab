function test_suite = PermutationCharacterTableTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_petersen
% Symmetric group of order 3
    ct = replab.ct.PermutationCharacterTable(replab.S(3));
    %vals = [1,1,1;2,0,-1;1,-1,1];
    %assert(isequal(ct.characterValues, vals))
end