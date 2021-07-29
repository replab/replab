function test_suite = PhasedMatrixPartitionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_chsh
    I = [1 2 3 4 5 % index matrix
         2 1 6 7 8
         3 6 1 9 10
         4 7 9 1 11
         5 8 10 11 1];
    % symmetry group
    g1 = [1 -2 -3 -4 -5];
    g2 = [1 4 5 2 3];
    g3 = [1 3 2 4 -5];
    pmp1 = replab.equi.PhasedMatrixPartition.fromIndexMatrix(I);
    G = replab.SignedPermutationGroup.of(g1,g2,g3);
    mu = replab.perm.GeneralizedSymmetricGroup.morphismFromSignedPermutationGroup(G);
    rep = G.naturalRep;
    h1 = mu.imageElement(g1); h2 = mu.imageElement(g2); h3 = mu.imageElement(g3);
    pmp2 = replab.equi.PhasedMatrixPartition.fromGeneralizedPermutations(2, {h1 h2 h3}, {h1 h2 h3});
    pmp3 = replab.equi.PhasedMatrixPartition.intersection(pmp1, pmp2);
    phase = [0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 1; 0 0 0 0 0; 0 0 1 0 0];
    block = [1 0 0 0 0; 0 1 0 2 2; 0 0 1 2 2; 0 2 2 1 0; 0 2 2 0 1];
    pmp4 = replab.equi.PhasedMatrixPartition.fromPhaseAndSubsetIndexMatrices(2, phase, block);
    assert(pmp3 == pmp4);
end
function test_intersection
    for i = 1:100
        pmp1 = replab.equi.PhasedMatrixPartition.fromIndexMatrix(randi([0 100],20,20));
        pmp2 = replab.equi.PhasedMatrixPartition.fromIndexMatrix(randi([0 100],20,20));
        pmp12 = replab.equi.PhasedMatrixPartition.intersection(pmp1, pmp2);
        pmp21 = replab.equi.PhasedMatrixPartition.intersection(pmp2, pmp1);
        assert(pmp12 == pmp21);
    end
end
function test_bug_marie
    subsetIndex = [0 2 1
                   1 3 0
                   2 0 4];
    phase = zeros(3, 3);
    pmp1 = replab.equi.PhasedMatrixPartition.fromIndexMatrix(subsetIndex);
    pmp2 = replab.equi.PhasedMatrixPartition.antisymmetric(3);
    pmp12 = replab.equi.PhasedMatrixPartition.intersection(pmp1, pmp2);
    subsetIndex_res = [0  1 -1
                       -1 0  0
                        1 0  0];
    pmp12_res = replab.equi.PhasedMatrixPartition.fromIndexMatrix(subsetIndex_res);
    assert(pmp12 == pmp12_res);
end