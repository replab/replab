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
    ctC = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 4);
end

function test_holt_Exercise2
    nGenerators = 2;
    relators = {[1 1 2 2], [-2 1 2 -1 -1 -1]};
    y = {};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 8);
end

function test_quaternion_presentation
    a = [2 6 1 8 7 3 4 5];
    b = [4 7 8 6 1 5 3 2];
    G = replab.PermutationGroup.of(a, b);
    relators = replab.fp.relatorsForPermutationGroup(G);
    assert(all(cellfun(@(r) G.isIdentity(G.imageLetters(r)), relators)));
end

function test_symmetric_group_presentation
    G = replab.SymmetricGroup(6);
    relators = replab.fp.relatorsForPermutationGroup(G);
    assert(all(cellfun(@(r) G.isIdentity(G.imageLetters(r)), relators)));
end

function test_dihedral_group_presentation
    if ReplabTestParameters.onlyFastTests
        return
    end
    G = replab.DihedralGroup(20);
    relators = replab.fp.relatorsForPermutationGroup(G);
    assert(all(cellfun(@(r) G.isIdentity(G.imageLetters(r)), relators)));
end

function test_Mathieu_M11
    if ReplabTestParameters.onlyFastTests
        return
    end
    r1 = [1 1];
    r2 = [2 2 2 2];
    r3 = [1 2 1 2 1 2 1 2 1 2  1 2 1 2 1 2 1 2 1 2  1 2];
    r4 = [1 2 2 1 2 2 1 2 2  1 2 2 1 2 2 1 2 2];
    r5 = [1 2 1 2 1 -2 1 2 1 2 2 1 -2 1 2 1 -2 1 -2];
    ctR = replab.fp.CosetTable.cosetEnumerationR(2, {r1 r2 r3 r4 r5}, {});
    ctC = replab.fp.CosetTable.cosetEnumerationC(2, {r1 r2 r3 r4 r5}, {});
    assert(isequal(ctR.C, ctC.C));
end

function test_Mathieu_M12_Holt
    if ReplabTestParameters.onlyFastTests
        return
    end
    nGenerators = 3;
    relators = {ones(1,11), [2,2], [3,3], [1,2,1,2,1,2], [1,3,1,3,1,3], ...
                repmat([2,3],1,10), [1,1,2,3,2,3,1,-3,-2,-3,-2]};
    y = {};
    ctR = replab.fp.CosetTable.cosetEnumerationR(nGenerators, relators, y);
    ctC = replab.fp.CosetTable.cosetEnumerationC(nGenerators, relators, y);
    assert(isequal(ctR.C, ctC.C));
    assert(size(ctR.C, 1) == 95040);
end

function test_Mathieu_M12_atlas
    if ReplabTestParameters.onlyFastTests
        return
    end
    G = replab.AbstractGroup.parsePresentation('< x, y | x^2 = y^3 = (x*y)^11 = [x,y]^6 = (x*y*x*y*x*y^-1)^6 = [x,y*x*y]^5 = 1 >');
    assert(G.order == 95040);
end