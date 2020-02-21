function [B_internal E_internal] = sageSpechtStandardBasis(d)
% Returns a basis for the standard representation embedded in the defining representation of S(d)
%
% Inspired by the Specht basis for irreducible representations from Sage Math
% (but not exactly the same).
%
% Args:
%   d (integer): Domain size
%
% Returns
% -------
% B_internal:
%   double(d,d-1), may be sparse: Basis
% E_internal:
%   double(d-1,d), may be sparse: Embedding map
    B_internal = ones(d,d-1)/d;
    for i = 1:d-1
        B_internal(i,i) = -(d-1)/d;
    end
    B_internal(:,1:2:d-1) = -B_internal(:,1:2:d-1);
    E_internal = sparse([1:(d-1) 1:(d-1)], [1:(d-1) ones(1,d-1)*d], [-ones(1,(d-1)) ones(1,(d-1))], d-1, d);
    E_internal(1:2:d-1,:) = -E_internal(1:2:d-1,:);
end
