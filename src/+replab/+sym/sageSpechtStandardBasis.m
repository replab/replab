function [injectiond projection] = sageSpechtStandardBasis(d)
% Returns a change of basis matrix for the standard representation embedded in the defining representation of S(d)
%
% Inspired by the Specht basis for irreducible representations from Sage Math
% (but not exactly the same).
%
% Args:
%   d (integer): Domain size
%
% Returns
% -------
% injectiond:
%   double(d,d-1), may be sparse: Injection map multiplied by d
% projection:
%   double(d-1,d), may be sparse: Embedding map
    injectiond = ones(d,d-1);
    for i = 1:d-1
        injectiond(i,i) = -(d-1);
    end
    injectiond(:,1:2:d-1) = -injectiond(:,1:2:d-1);
    projection = sparse([1:(d-1) 1:(d-1)], [1:(d-1) ones(1,d-1)*d], [-ones(1,(d-1)) ones(1,(d-1))], d-1, d);
    projection(1:2:d-1,:) = -projection(1:2:d-1,:);
end
