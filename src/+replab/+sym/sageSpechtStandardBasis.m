function [H_internal F_internal] = sageSpechtStandardBasis(d)
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
% H_internal:
%   double(d,d-1), may be sparse: Basis
% F_internal:
%   double(d-1,d), may be sparse: Embedding map
    H_internal = ones(d,d-1)/d;
    for i = 1:d-1
        H_internal(i,i) = -(d-1)/d;
    end
    H_internal(:,1:2:d-1) = -H_internal(:,1:2:d-1);
    F_internal = sparse([1:(d-1) 1:(d-1)], [1:(d-1) ones(1,d-1)*d], [-ones(1,(d-1)) ones(1,(d-1))], d-1, d);
    F_internal(1:2:d-1,:) = -F_internal(1:2:d-1,:);
end
