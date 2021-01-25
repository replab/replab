function [G, images] = monomialImages(rep, preimages)
% Returns the phase order when the representation is monomial, or ``0`` when it cannot be shown to be monomial
%
% Args:
%   rep (`+replab.Rep`): Representation to extract the images of
%
% Returns
% -------
%   G: `+replab.perm.GeneralizedSymmetricGroup`
%     Group of monomial images, or ``[]`` if the representation is not monomial
%   images: cell(1,\*) of elements of ``G``
%     Requested images
    if isa(rep, 'replab.rep.DerivedRep')
        [G, images] = replab.rep.monomialImages(rep.parent, preimages);
        if isempty(G)
            return
        end
        if xor(rep.conjugate, rep.transpose)
            m = G.m;
            % apply conjugate
            images = cellfun(@(g) [g(1,:); m+1-g(2,:)], images, 'uniform', 0);
        end
    elseif isa(rep, 'replab.rep.TrivialRep')
        n = rep.dimension;
        G = replab.perm.GeneralizedSymmetricGroup(n, 1);
        images = repmat({[1:n; zeros(1, n)]}, 1, length(preimages));
    elseif isa(rep, 'replab.rep.RepByImages_monomial')
        G = rep.morphism.target;
        images = cellfun(@(x) rep.morphism.imageElement(x), preimages, 'uniform', 0);
    else
        G = [];
        images = [];
    end
end
