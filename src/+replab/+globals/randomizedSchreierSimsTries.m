function value = randomizedSchreierSimsTries(newValue)
% Gets/sets the number of sifted elements before the BSGS chain is declared complete
%
% This is the number of successive failed attempts to generate a new strong generator
% before deciding the chain is complete in the randomized Schreier-Sims algorithm;
% the probability of failure is then less than 1/2^value.
%
% Args:
%   newValue (integer, optional): New value
%
% Returns:
%   integer: The current value.
    persistent RandomizedSchreierSimsTries;
    if nargin == 1
        RandomizedSchreierSimsTries = newValue;
    elseif isempty(RandomizedSchreierSimsTries)
        RandomizedSchreierSimsTries = 100;
    end
    value = RandomizedSchreierSimsTries;
end
