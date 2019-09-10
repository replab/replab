function shouldProduceAnError(func, nbOutputs)
% shouldProduceAnError(func, [nbOutputs])
%
% Checks that the action of function 'func' produces an error. In this case
% the error is silenced. If the function doesn't produce an error, an error
% is thrown.
%
% nbOutputs is the number of desired outcomes for the function (0 by
% default)

if nargin < 1
    error('not enough arguments');
end

if nargin < 2
    nbOutputs = 0;
end

switch nbOutputs
    case 0
        try
            func();
            error('The error test failed')
        catch me
            if isequal(me.message, 'The error test failed')
                assert(false);
            end
        end
    case 1
        try
            [a] = func();
            error('The error test failed')
        catch me
            if isequal(me.message, 'The error test failed')
                assert(false);
            end
        end
    case 2
        try
            [a b] = func();
            error('The error test failed')
        catch me
            if isequal(me.message, 'The error test failed')
                assert(false);
            end
        end
    case 3
        try
            [a b c] = func();
            error('The error test failed')
        catch me
            if isequal(me.message, 'The error test failed')
                assert(false);
            end
        end
    otherwise
        error('unsupported call to shouldProduceAnError');
end

end
