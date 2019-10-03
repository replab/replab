function test_suite = RepresentationsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    
    % test orbit decomposition
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.definingRep;
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
    % test standard representation
    G = replab.Permutations(5);
    Gstd = G.standardRep;
    test_suite = replab.RepLaws(Gstd).addTestCases(test_suite);
    % test quaternion type representations
    Q = replab.signed.Permutations.quaternionGroup;
    S3 = replab.S(3);
    W = S3.wreathProduct(Q);
    rho = W.primitiveRep(Q.definingRep);
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
    rho = W.imprimitiveRep(Q.definingRep);
    I = rho.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
end

function test_symmetric_group_representations
    G = replab.Permutations(4).subgroup({[2 1 3 4] [2 3 4 1]});
    rho = G.definingRep;
    I = rho.decomposition;
    assertEqual(I.nComponents, 2);
    assertEqual(cellfun(@(c) c.irrepDimension, I.components), [1 3]);
    assertEqual(cellfun(@(c) c.multiplicity, I.components), [1 1]);
end

function test_representation_of_cyclic_group
    C12 = replab.Permutations(12).cyclicSubgroup;
    rep = C12.definingRep;
    dec = rep.decomposition;
    d = cellfun(@(iso) iso.irrepDimension, dec.components);
    m = cellfun(@(iso) iso.multiplicity, dec.components);
    t = cellfun(@(iso) iso.irrep(1).irrepInfo.divisionAlgebra, dec.components);
    assertEqual(d, [1 1 2 2 2 2 2]);
    assertEqual(m, [1 1 1 1 1 1 1]);
    assertEqual(t, 'RRCCCCC');
end
