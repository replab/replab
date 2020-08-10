function test_suite = graphAutomorphisms()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

end

function test_Frucht
    MFrucht = [0  1  0  0  0  0  0  0  0  0  0  1
               1  0  1  1  0  0  0  0  0  0  0  0
               0  1  0  0  0  0  0  0  0  0  1  0
               0  1  0  0  1  1  0  0  0  0  0  0
               0  0  0  1  0  1  0  1  0  0  0  0
               0  0  0  1  1  0  0  0  0  1  0  0
               0  0  0  0  0  0  0  1  1  0  0  1
               0  0  0  0  1  0  1  0  1  0  0  0
               0  0  0  0  0  0  1  1  0  1  0  0
               0  0  0  0  0  1  0  0  1  0  1  0
               0  0  1  0  0  0  0  0  0  1  0  1
               1  0  0  0  0  0  1  0  0  0  1  0];
    graph = replab.UndirectedGraph.fromAdjacencyMatrix(MFrucht);
    automorphismGroup = graph.automorphismGroup;
    assert(automorphismGroup.order == 1);
end

function test_Coxeter
    % The Coxeter graph
    MCoxeter = [0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0
                1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
                0  1  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  1  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
                0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0
                0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1
                0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0
                0  0  0  1  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1  0
                0  0  1  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  1
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  1  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  1  0  0
                0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0
                0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  1
                0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0
                0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  1  0  0  0
                1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  0
                0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0
                1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
                0  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0
                0  0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0];
    graph = replab.UndirectedGraph.fromAdjacencyMatrix(MCoxeter);
    automorphismGroup = graph.automorphismGroup;
    assert(automorphismGroup.order == 336);
end

function test_dodecahedron
    % The dodecahedron admits 120 automorphisms
    MDodecahedron = [0  1  0  0  0  0  0  0  0  1  1  0  0  0  0  0  0  0  0  0
                     1  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
                     0  1  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0
                     0  0  1  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0
                     0  0  0  1  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0
                     0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0
                     0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  1  0  0  0
                     0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  1  0  0
                     0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  0  0  0  1  0
                     1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1
                     1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1  0
                     0  1  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  1
                     0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0  0
                     0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0  0
                     0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0  0
                     0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0  0
                     0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1  0
                     0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  1
                     0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  1  0  0  0
                     0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0  1  0  0];
    graph = replab.UndirectedGraph.fromAdjacencyMatrix(MDodecahedron);
    automorphismGroup = graph.automorphismGroup;
    assert(automorphismGroup.order == 120);
end

