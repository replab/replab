function value = doubleSdpTol(newValue)
    persistent DoubleSdpTol;
    if nargin == 1
        DoubleSdpTol = newValue;
    elseif isempty(DoubleSdpTol)
        DoubleSdpTol = 1e-5;
    end
    value = DoubleSdpTol;
end
