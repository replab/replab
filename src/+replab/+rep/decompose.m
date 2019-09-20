function irreps = decompose(rep)
% Decomposes the given representation into irreducible subrepresenations
%
% Note: see the default methods applied in `replab.dispatchDefaults`
%
% Args:
%   rep (replab.SubRep): (Part of) the representation to decompose
%
% Returns:
%   row cell array of replab.SubRep: A cell array of irreducible subrepresentations
    assert(isa(rep, 'replab.SubRep'));
    irreps = replab.dispatch('call', 'replab.rep.decompose', rep);
end
