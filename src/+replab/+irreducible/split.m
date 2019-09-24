function sub = split(rep, samples, U)
% Decomposes the given representation into subrepresentations
%
% Note: see the default methods applied in `replab.dispatchDefaults`
%
% Args:
%   rep (replab.Rep): Representation to decompose
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples
%   U (double matrix, can be sparse): Orthonormal row vectors describing a subrepresentation of `rep`
%
% Returns:
%   row cell array of replab.SubRep: A cell array of irreducible subrepresentations
%                                    with trivial representations identified with the label '1'
    assert(isa(rep, 'replab.Rep'));
    switch nargin
      case 2
        replab.irreducible.tell('dispatch split dimension');
        sub = replab.dispatch('call', 'replab.irreducible.split', rep, samples);
      case 3
        replab.irreducible.tell('dispatch split dimension %d', size(U, 1));
        sub = replab.dispatch('call', 'replab.irreducible.split', rep, samples, U);
    end
end
