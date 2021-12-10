function res = multikron(factors, type)
% Computes the Kronecker product of multiple factors with a given type
%
% Args:
%   factors (cell(1,\*) of double(\*,\*) or `+replab.cyclotomic` (\*,\*)): Factors
%   type ('double' 'double/sparse' or 'exact', optional): Type of the result, default: double
    if isempty(factors)
        if strcmp(type, 'exact')
            res = replab.cyclotomic.eye(1);
        else
            res = 1;
        end
    else
        res = factors{1};
        for i = 2:length(factors)
            res = kron(res, factors{i});
        end
        if strcmp(type, 'exact')
            assert(isa(res, 'replab.cyclotomic'));
        else
            res = double(res);
        end
        if strcmp(type, 'double')
            res = full(res);
        end
    end
end