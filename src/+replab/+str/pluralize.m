function s = pluralize(n, singular, plural)
% Returns a string description of a given number of given things
%
% ``pluralize(3, 'horse', 'horses')`` is ``'three horses'``
% ``pluralize(0, 'apple', 'apples')`` is ``'zero apple'`` and so on.
%
% The plural form can be omitted if it is formed by the suffix -s.
    if nargin < 3
        plural = [singular 's'];
    end
    numbers = {'one' 'two' 'three' 'four' 'five' 'six' 'seven' 'eight' 'nine' 'ten'};
    if n == 0
        s = sprintf('zero %s', singular);
    elseif n == 1
        s = sprintf('one %s', singular);
    elseif n <= 10
        s = [numbers{n} ' ' plural];
    else
        s = [num2str(n) ' ' plural];
    end
end
