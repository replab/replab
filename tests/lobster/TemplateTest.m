function test_suite = TemplateTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_empty(self)
    tpl = replab.lobster.Template('');
    assertEqual(tpl.render(), '');
end

function test_simple(self)
    tpl = replab.lobster.Template('This is a test string.');
    assertEqual(tpl.render(), 'This is a test string.');
end

function test_int_var(self)
    context.var = 1;
    tpl = replab.lobster.Template('{{ var }}');
    assertEqual(tpl.render(context), '1');
end

function test_string_var(self)
    context.var = 'stuff';
    tpl = replab.lobster.Template('{{ var }}');
    assertEqual(tpl.render(context), 'stuff');
end

function test_text_and_var(self)
    context.var = 1;
    tpl = replab.lobster.Template('This is {{ var }}');
    assertEqual(tpl.render(context), 'This is 1');
end

function test_var_and_text(self)
    context.var = 1;
    tpl = replab.lobster.Template('{{ var }} is cool');
    assertEqual(tpl.render(context), '1 is cool');
end

function test_var_with_map_access(self)
    if ~replab.compat.isOctave % old versions of Octave don't have containers
        context.var = containers.Map('some_key', 'the value');
        tpl = replab.lobster.Template('{{ var(''some_key'') }} is cool');
        assertEqual(tpl.render(context), 'the value is cool');
    end
end

function test_undefined_var_error(self)
    context = struct();
    tpl = replab.lobster.Template('{{ var }} is cool');
    if ~replab.compat.isOctave % old versions of Octave rather throw Octave:invalid-fun-call
        assertExceptionThrown(@() tpl.render(), 'Lobster:TemplateContextError');
    end
end

function test_if_true_with_no_context(self)
    tpl = replab.lobster.Template('{% if true %}You should see this{% endif %}');
    assertEqual(tpl.render(), 'You should see this');
end

function test_if_false_with_no_context(self)
    tpl = replab.lobster.Template('{% if false %}You should not see this{% end %}');
    assertEqual(tpl.render(), '');
end

function test_if_true_with_else(self)
    tpl = replab.lobster.Template('{% if true %}Show this{% else %}Not this{% end %}');
    assertEqual(tpl.render(), 'Show this');
end

function test_if_false_with_else(self)
    tpl = replab.lobster.Template('{% if false %}Show this{% else %}Not this{% end %}');
    assertEqual(tpl.render(), 'Not this');
end

function test_if_with_conditional(self)
    tpl = replab.lobster.Template('{% if length(1:5) > 4 %}You should see this{% endif %}');
    assertEqual(tpl.render(), 'You should see this');
end

function test_for_with_array(self)
    tpl = replab.lobster.Template('{% for k in 1:5 %}{{ k }} {% end %}');
    assertEqual(tpl.render(), '1 2 3 4 5 ');
end

function test_for_with_empty_array(self)
    tpl = replab.lobster.Template('{% for k in [] %}{{ k }} {% end %}');
    assertEqual(tpl.render(), '');
end

function test_for_with_cell(self)
    context.collection = num2cell(1:5);
    tpl = replab.lobster.Template('{% for k in 1:5 %}{{ k }} {% end %}');
    assertEqual(tpl.render(context), '1 2 3 4 5 ');
end

function test_for_with_empty_cell(self)
    context.collection = cell(0);
    tpl = replab.lobster.Template('{% for k in collection %}{{ k }} {% end %}');
    assertEqual(tpl.render(context), '');
end

function test_for_with_struct(self)
    context.collection = struct('val', num2cell(1:5));
    tpl = replab.lobster.Template('{% for k in collection %}{{ k.val }} {% end %}');
    assertEqual(tpl.render(context), '1 2 3 4 5 ');
end

function test_for_with_empty_struct(self)
    context.collection = struct([]);
    tpl = replab.lobster.Template('{% for k in collection %}{{ k.val }} {% end %}');
    assertEqual(tpl.render(context), '');
end
