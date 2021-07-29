function mu = relabelings(n, m, k, convention)
% Returns the permutation group that relabels a probability distribution/Bell inequality
%
% This function provides the permutations using two labeling conventions:
%
% * The ``Pabxy`` convention assumes that the permutation is applied to a multi-dimensional
%   array whose coefficients are accessed as ``P_AB(a,b,x,y)`` or ``P_AB(a,b,c,x,y,z)``, etc, where
%   ``a``, ``b``, ``c`` are output indices and ``x``, ``y``, ``z`` are input indices.
%   Of course, we cannot apply a relabeling to a multi-dimensional array, so the permutation
%   is made on the flattened vector ``P_AB(:)``.
%
% * The ``kron`` convention assumes that the vector representing ``P_AB`` obeys that
%   ``P_AB = kron(P_A, P_B)``, and that, for the vector ``P_A``, the coefficients correspond
%   to the indices ``(a,x)`` enumerated in the order ``(1,1), (2,1), ..., (d,1), (1,2), ..., (d,m)``.
%   Note that a vector in the ``kron`` convention corresponds, after reshape, to a multi-dimensional
%   array with indices ``P_AB(b,y,a,x)``.
%
% Args:
%   n (integer): Number of parties
%   m (integer): Number of inputs
%   k (integer): Number of outputs
%   convention ('Pabxy', 'kron', optional): Labeling convention, default: ``'Pabxy'``
%
% Returns:
%   `+replab.FiniteIsomorphism`: Isomorphism from a structured to a permutation group that represents all relabelings
    if nargin < 4 || isempty(convention)
        convention = 'Pabxy';
    end
    ind = 1:(k*m)^n;
    switch convention
      case 'Pabxy'
        ind = reshape(1:(k*m)^n, [ones(1,n)*k ones(1,n)*m]);
        order = [n:-1:1
                 n+(n:-1:1)];
        ind = permute(ind, order(:)');
        ind = ind(:)';
      case 'kron'
        % do nothing
      otherwise
        error('Invalid convention %s', convention);
    end
    d = length(ind);
    Sd = replab.S(d);
    generators = cell(1, 0);
    Sparties = replab.S(n);
    Sinputs = replab.S(m);
    Soutputs = replab.S(k);
    W1 = Sinputs.wreathProduct(Soutputs); % single party
    W = Sparties.wreathProduct(W1);
    mu = W.isomorphismByFunction(replab.S(d), @(w) Sd.leftConjugate(ind, W.primitivePermutation(w, @(w1) W1.imprimitivePermutation(w1))));
end
