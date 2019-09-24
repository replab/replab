function sub = split(rep, samples, U)
% Decomposes the given representation into subrepresentations
%
% Note: see the default methods applied in `replab.dispatchDefaults`
%
% The dispatch function tests whether `U` is the identity and discards the argument if that's the case.
%
% Args:
%   rep (replab.Rep): Representation to decompose
%   samples (replab.rep.OnDemandSamples, optional): Lazy evaluation of various samples
%   U (double matrix, can be sparse, optional): Orthonormal row vectors describing a subrepresentation of `rep`
%
% Returns:
%   row cell array of replab.SubRep: A cell array of irreducible subrepresentations
    assert(isa(rep, 'replab.Rep'));
    if nargin == 1 || isempty(samples)
        samples = replab.rep.OnDemandSamples(rep);
    end
    if nargin < 3 || isempty(U) || replab.iseye(U)
        sub = replab.dispatch('call', 'replab.rep.split', rep, samples);
    else
        sub = replab.dispatch('call', 'replab.rep.split', rep, samples, U);
    end
end
