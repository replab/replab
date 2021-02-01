function value = bsgsFailureProbability(newValue)
    persistent BsgsFailureProbability;
    if nargin == 1
        BsgsFailureProbability = newValue;
    elseif isempty(BsgsFailureProbability)
        BsgsFailureProbability = 2^-100;
    end
    value = BsgsFailureProbability;
end
