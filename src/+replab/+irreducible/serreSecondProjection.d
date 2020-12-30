function P = serreSecondProjection(rep, irrep, alpha, beta)
% Computes Serre's projector P_{\alpha\beta} for a given explicit form of an irreducible representation
%
% For more background see:
%
% * Serre's textbook, Chapter 2
%
% * Faessler and Stiefel Chapter 5 on how to generate symmetry adapted basis, in which a
%   matrix with the symmetry of the representation becomes block-diagonal.
%
% The representing matrices for the ``irrep`` can be provided either by a
% RepLAB decomposition or from the literature, if the irreps of a given
% representation are known.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%   irrep (`+replab.Rep`): Irreducible representation
%
% Returns:
%   double(\*,\*), may be sparse: Projector of dimension ``rep.dimension x rep.dimension``
%
    id = irrep.dimension;
    rd = rep.dimension;
    P = cell(id, id);
    G = rep.group;
    elements = G.elements.toCell;

    s = zeros(rd, rd);
    for i = 1:length(elements)
        t = elements{i};
        repimg = irrep.image(G.inverse(t)); % Returns the image in the irrep of the inverse of an element t in the group G
        rep_ab = repimg(alpha, beta); % Follow the notation from Faessler and Stiefel on how to construct the projectors
        s = s + rep_ab * rep.image(t); % calculating the sum in the projection formula (Serre, Chapter 2.7); for this we need the representing matrix of the element in the full representation, thus making the linear map of dimensions equal to the full representation.
    end
    P = s * id / double(G.order); % returns the projection map
end
