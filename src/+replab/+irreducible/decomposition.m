function I = decomposition(rep)
% Decomposes the given representation into an irreducible decomposition, identifying details
%
% Details include grouping copies (multiplicities) in isotypic components, identifying division
% algebras (for representations over the reals), and making sure that complex-type and quaternion-type
% real representations are expressed in a canonical basis for the division algebra encoding.
%
% Note: see the default methods applied in `+replab.dispatchDefaults`
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   `+replab.Irreducible`: The irreducible decomposition
    replab.irreducible.tell('dispatch decomposition')
    assert(isa(rep, 'replab.Rep'));
    I = replab.dispatch('call', 'replab.irreducible.decomposition', rep);
end
