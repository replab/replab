function value = doubleEigTol(newValue)
    persistent DoubleEigTol;
    if nargin == 1
        DoubleEigTol = newValue;
    elseif isempty(DoubleEigTol)
        DoubleEigTol = 1e-10;
    end
    value = DoubleEigTol;
end
