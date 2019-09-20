function sub = decompose(rep)
% Decomposes the given representation into subrepresenations (one step)
%
% Note: see the default methods applied in `replab.dispatchDefaults`
%
% This method needs to be called recursively to complete the process
%
% Args:
%   rep (replab.SubRep): (Part of) the representation to decompose
%
% Returns:
%   row cell array of replab.SubRep: A cell array of subrepresentations, not necessarily irreducible
    assert(isa(rep, 'replab.SubRep'));
    sub = replab.dispatch('call', 'replab.rep.decompose', rep);
end
