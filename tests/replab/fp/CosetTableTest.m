function test_suite = CosetTableTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_Holt_5_4
    relators = {[1 2 -3], [2 3 -4], [3 4 -5], [4 5 -1], [5 1 -2]};
    ctR = replab.fp.CosetTable.cosetEnumerationR(5, relators, {[1]});
    ctC = replab.fp.CosetTable.cosetEnumerationC(5, relators, {[1]});
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 1);
    relators = {[1 2 -3], [2 3 -4], [3 4 -5], [4 5 -1], [5 1 -2]};
    ctR = replab.fp.CosetTable.cosetEnumerationR(5, relators, {});
    ctC = replab.fp.CosetTable.cosetEnumerationC(5, relators, {});
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 11);
end

function test_Holt_5_1
    nGenerators = 2;
    relators = {[1 1 1], [2 2 2], [-1 -2 1 2]};
    y = {[1]};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 3);
end

function test_Holt_5_2
    nGenerators = 2;
    relators = {[1 1], [2 2 2], [1 2 1 2 1 2]};
    y = {[1 2]};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 4);
end

function test_Holt_5_3
    nGenerators = 2;
    relators = {[1 1 2 2], [1 1 1 2 2 2 2 2]};
    y = {};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 4);
end

function test_holt_Exercise2
    nGenerators = 2;
    relators = {[1 1 2 2], [-2 1 2 -1 -1 -1]};
    y = {};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 8);
end

function test_holt_M12
    nGenerators = 3;
    relators = {ones(1,11), [2,2], [3,3], [1,2,1,2,1,2], [1,3,1,3,1,3], ...
                repmat([2,3],1,10), [1,1,2,3,2,3,1,-2,-3,-2,-3]};
    y = {};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 1);
end