function test_suite = helpTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    % Note: here we mostly generate some code coverage
end

function test_help_external
    assert(~isempty(evalc('help plot')));
end

function test_help_package
    assert(~isempty(evalc('help replab')));
end

function test_help_class
    assert(~isempty(evalc('help replab.Monoid')));
    assert(~isempty(evalc('help -f replab.Monoid')));
    assert(~isempty(evalc('help replab.Monoid.identity')));
    assert(~isempty(evalc('help -f replab.Monoid.identity')));
    assert(~isempty(evalc('help replab.Monoid.compose')));
    assert(~isempty(evalc('help -f replab.Monoid.compose')));
end

function test_help_function
    assert(~isempty(evalc('help replab.S')));
    assert(~isempty(evalc('help -f replab.S')));
end
