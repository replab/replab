function test_suite = sizeTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_general
    global matrix231 matrix23451 matrix23451H
    matrix = matrix23451;
    s12 = size(matrix);
    assert(s12(1) == size(matrix,1));
    assert(s12(2) == size(matrix,2));
    [s1 s2] = size(matrix);
    assert(s1 == size(matrix,1));
    assert(s2 == size(matrix,2));
end

function test_inputs
    global matrix231 matrix23451 matrix23451H
    matrix = matrix231;
    
    % Octave > 4.2 has some trouble with anonymous functions that involve class objects
    isOctave = replab.settings.isOctave;
    ver = version;
    
    if ~isOctave || (isOctave && isequal(ver(1:3),'4.2'))
        shouldProduceAnError(@(x) size(matrix, 3));
    end
    
    shouldProduceAnError(@(x) size(1, matrix));
end
