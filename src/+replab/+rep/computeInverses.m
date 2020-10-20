function inverseImages = computeInverses(source, targetCompose, preimages, images)
% Given a list of preimages and images defining a homomorphism, efficiently computes the inverse images
%
% This function avoids computing inverse images directly when it is time consuming to do so. Rather,
% it identifies existing inverse relationships in the preimages to identify inverse images directly,
% and uses the preimages element order to compute the inverse image by taking a power of it.
%
% Args:
%   source (`+replab.FiniteGroup`): Source group where inverse computations are fast
%   targetCompose (function_handle): Function taking two arguments that computes composition in the target group
%   preimages (cell(1,\*) of elements of ``source``): Preimages
%   images (cell(1,\*) of elements of the target group): Images
%
% Returns:
%   cell(1,\*) of elements of the target group: Inverse images
    n = length(preimages);
    assert(length(images) == n);
    invInd = replab.mrp.inverseIndices(source, preimages);
    inverseImages = cell(1, n);
    for i = 1:n
        if invInd(i) ~= 0
            inverseImages{i} = images{invInd(i)};
        else
            o = source.elementOrder(preimages{i});
            inverseImages{i} = replab.util.repeatedSquaring(images{i}, o-1, targetCompose);
        end
    end
end
