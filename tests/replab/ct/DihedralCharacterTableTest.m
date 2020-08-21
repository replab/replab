function test_suite = DihedralCharacterTableTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_D3
% Dihedral group of order 6
    ct = replab.ct.DihedralCharacterTable(3);
    %    expressions = {'1','1','1';'1','-1','1';'2','0','E(3)^1+E(3)^-1'};
    %assert(isequal(ct.characterExpressions, expressions))
end

function test_D4
% Dihedral group of order 8
    ct = replab.ct.DihedralCharacterTable(4);
    %expressions = {'1','1','1','1','1';
    %               '1','-1','-1','1','1';
    %               '1','-1','1','-1','1';
    %               '1','1','-1','-1','1';
    %               '2','0','0','E(4)^1+E(4)^-1','E(4)^2+E(4)^-2'};
    %assert(isequal(ct.characterExpressions, expressions))
end