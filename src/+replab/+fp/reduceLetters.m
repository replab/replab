function x = reduceLetters(x)
% Finds the reduced form of the word described the given letters
%
% It eliminates all products of the form ``[i -i]`` or ``[-i i]`` for the index ``i`` of a generator.
%
% Args:
%   x (integer(1,\*)): Letters
%
% Returns:
%   integer(1,\*): Letters of the reduced word
    if isempty(x)
        x = zeros(1, 0);
        return
    end
    i = find(x(1:end-1) == -x(2:end));
    if ~isempty(i)
        while i < length(x)
            if x(i) == -x(i+1)
                x = [x(1:i-1) x(i+2:end)];
                if i > 1
                    i = i - 1;
                end
            else
                i = i + 1;
            end
        end
    end
end
