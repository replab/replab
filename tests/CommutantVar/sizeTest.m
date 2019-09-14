function test_suite = sizeTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_oneGroup
    % We do a sanity check with one group
    global matrix231 matrix23451
    matrix = matrix23451;
    s12 = size(matrix);
    assert(s12(1) == size(matrix,1));
    assert(s12(2) == size(matrix,2));
end

function test_inputs
    global matrix231 matrix23451
    matrix = matrix231;
    
    % Octave > 4.2 has some trouble with anonymous functions that involve class objects
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    ver = version;
    
    if ~isOctave || (isOctave && isequal(ver(1:3),'4.2'))
        shouldProduceAnError(@(x) size(matrix, 3));
    end
    
    shouldProduceAnError(@(x) size(1, matrix));
end
