function [orbit, transversal] = orbitTransversal(degree, generators, alpha)
% Returns the transversal for the orbit of a point
%
% A transversal for the orbit ```O = { g(alpha) : g \in G }`` is given by a set
% ``{ g_beta }`` for ``g_beta(alpha) = beta`` and ``beta \in O``
%
% Args:
%   generators (integer(\*,\*)): Group generators given as matrix columns
%   b (integer): Point to compute the orbit of
%   ordering (integer(1,\*), optional): Ordering, default value ``[1 2 ... n]``
%
% Returns
% -------
% orbit:
%   integer(1,\*): Orbit of ``alpha``
% transversal:
%   integer(\*,\*): Transversal elements as columns, in the order of the returned ``orbit``
    nG = size(generators, 2);
    orbit = [alpha];
    % see the algorithm in `+replab.bsgs.Chain.completeOrbit`
    iOrbit = zeros(1, degree);
    iOrbit(alpha) = 1;
    transversal = [1:degree]';
    toTest = orbit;
    while ~isempty(toTest)
        imgs = generators(toTest, :);
        mask = iOrbit(imgs) == 0;
        mask = reshape(mask, size(imgs)); % to be sure
        [bInd, genInd] = find(mask);
        for j = 1:length(bInd)
            gen = generators(:, genInd(j));
            b = toTest(bInd(j)); % orbit element
            newb = gen(b);
            if iOrbit(newb) == 0
                ind = iOrbit(b);
                newt = gen(transversal(:,ind));
                k = length(orbit) + 1;
                orbit(1,k) = newb;
                iOrbit(newb) = k;
                transversal(:,k) = newt;
            end
        end
        toTest = imgs(mask);
        toTest = toTest(:)';
        if length(toTest) > 1
            toTest = sort(toTest);
            mask1 = [true, toTest(1:end-1) ~= toTest(2:end)];
            toTest = toTest(mask1);
        end
    end
end
