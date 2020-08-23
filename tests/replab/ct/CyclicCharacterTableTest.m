function test_suite = CyclicCharacterTableTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_C3
% Cyclic group or order 3
    ct = replab.ct.CyclicCharacterTable(3);
    %expressions = {'E(3)^0', 'E(3)^0', 'E(3)^0';
    %               'E(3)^0', 'E(3)^1', 'E(3)^2';
    %               'E(3)^0', 'E(3)^2', 'E(3)^1'};
    %assert(isequal(ct.characterExpressions, expressions))
end