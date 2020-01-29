function out = eval_with_context(expression, context)
% Evaluates an expression in context (?)
%
% TODO: reraise the error if necessary

    % Unpack the context stuct into this workspace
    fields__ = fieldnames(context);
    for k__ = 1:length(fields__)
        feval(@() assignin('caller', fields__{k__}, context.(fields__{k__})));
    end

    clear fields__ k__;

    try
        out = eval(expression);
    catch ME
        err_msg = {
            ''
            'The following error occurred whilst evaluating a matlab expression. '
            ''
            '   [%s] %s'
            ''
            'The expression was: '
            '    >> %s'
            ''
            'The context variable contained the following fields: '
            '    %s'
            ''
        };

        fields = strjoin(fieldnames(context), ', ');
        error('Lobster:TemplateContextError', strjoin(err_msg, '\n'), ...
            ME.identifier, ME.message, expression, fields);
    end


end
