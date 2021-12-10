function [P, n] = isotypicProjector(rep, varargin)
% Returns, for a given representation, the projector into the isotypic subspace corresponding to the given irrep
%
% This uses the projection formula from Theorem 8 of
% J.-P. Serre, Linear Representations of Finite Groups. Springer, 1977
%
% Note that the formula works as well when ``rep`` is a real representation.
%
% * For real-type real representations, it works as it is.
% * For complex-type real representations, the character is the sum of the conjugate complex representations, so the
%   projector will be the sum of the projectors into both conjugates.
% * For quaternion-type real representations, the character is the sum of two complex equivalent representations, so
%   the projector needs to be rescaled.
%
% Adapted from the code by Sotiris Mygdalas
%
% Args:
%   rep (`+replab.Rep`): Representation of dimension ``D``
%
% Keyword Args:
%   type ('exact', 'double', 'double/sparse'): Type of the projector to compute
%   irrep (`+replab.Rep`): Irreducible representation
%   irrepCharacter (`+replab.Character`): Irreducible character
%
% Returns
% -------
%   P: double(D,D) or `+replab.cyclotomic` (D,D), may be sparse
%     Projector
%   n: integer
%     Returns the multiplicity of the given irreducible character
    args = struct('irrep', [], 'irrepCharacter', [], 'type', 'double');
    args = replab.util.populateStruct(args, varargin);
    assert(isa(rep.group, 'replab.FiniteGroup'), 'The Serre projection formula only works for finite groups.');
    type = args.type;
    if strcmp(type, 'exact')
        assert(rep.isExact);
    end
    if isempty(args.irrepCharacter)
        assert(~isempty(args.irrep), 'One of the irrep or irrepCharacter keyword arguments must be specified');
        irrep = args.irrep;
        assert(rep.field == irrep.field);
        assert(irrep.isIrreducible);
        assert(rep.group == irrep.group);
        if strcmp(args.type, 'exact')
            irrepCharacter = replab.Character.fromRep(args.irrep);
            charFun = @(g) irrepCharacter.value(g);
        else
            charFun = @(g) full(trace(irrep.image(g, 'double/sparse')));
        end
    else
        assert(isempty(args.irrep), 'Only one of irrep or irrepCharacter can be specified');
        irrepCharacter = args.irrepCharacter;
        charFun = @(g) irrepCharacter.value(g);
    end
    D = rep.dimension;
    elements = rep.group.elements.toCell;
    P = [];
    for i = 1:length(elements)
        t = elements{i};
        if strcmp(type, 'exact')
            X = conj(charFun(t))*rep.image(t, 'exact');
        else
            X = conj(double(charFun(t)))*rep.image(t, 'double/sparse');
        end
        if i == 1
            P = X;
        else
            P = P + X;
        end
    end
    if strcmp(type, 'exact')
        P = P / replab.cyclotomic.fromVPIs(rep.group.order);
    else
        P = P / double(rep.group.order);
    end
    n = trace(P)/trace(P*P);
    P = P * n;
end
