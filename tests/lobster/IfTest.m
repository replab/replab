function test_suite = IfTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_if_false_with_context(self)
    values = {false, 0, '', []};
    for i = 1:length(values)
        falsy_value = values{i};
        context.var = falsy_value;
        tpl = replab.lobster.Template('{% if var %}You should not see this{% end %}');
        assertEqual(tpl.render(context), '');
    end
end

function test_if_true_with_context(self)
    values = {true, 1, -1, 2, 5, -7, 'true', 'stuff', [1, 1]};
    for i = 1:length(values)
        truthy_value = values{i};
        context.var = truthy_value;
        tpl = replab.lobster.Template('{% if var %}You should see this{% endif %}');
        assertEqual(tpl.render(context), 'You should see this');
    end
end
