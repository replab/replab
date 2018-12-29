function test_suite = RepresentationsTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_symmetric_group_representations
    G = replab.PermutationGroup.fromGenerators({[2 1 3 4] [2 3 4 1]});
    rho = G.naturalRepresentation;
    I = rho.irreducible;
    assertEqual(I.nComponents, 2);
    assertEqual(I.dimensions, [1 3]);
    assertEqual(I.multiplicities, [1 1]);
end
