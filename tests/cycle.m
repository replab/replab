function [input] = cycle(input, varargin)

%convention used - cycle of (123) corresponds to 1 -> 3, 2 -> 1, 3 -> 2

for j = 1:length(varargin)
    
    cycle = varargin{j};
    c = length(cycle);
    in = input(cycle(1));
    for i = 1:c-1
       input(cycle(i)) = input(cycle(i + 1));
    end
    input(cycle(c)) = in;
end

