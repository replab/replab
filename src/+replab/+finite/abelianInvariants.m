function a = abelianInvariants(group)
% See `.abelianInvariants`
%
% Inspired by GAP System "AbelianInvariants"
    if group.isTrivial
        a = zeros(1, 0);
        return
    end
    a = zeros(1, 0);
    G = group;
    dsg = group.derivedSubgroup;
    primeDivisors = unique(double(factor(group.order)));
    for i = 1:length(primeDivisors)
        p = primeDivisors(i);
        ranks = zeros(1, 0);
        r = 0;
        while r ~= 1
            H = dsg;
            for j = 1:G.nGenerators
                H = H.closure(G.composeN(G.generator(j), p));
            end
            r = G.order / H.order;
            G = H;
            if r ~= 1
                ranks = [ranks length(factor(r))];
            end
        end
        if ~isempty(ranks)
            l = ones(1, ranks(1));
            for i = ranks
                l(1:i) = l(1:i) * p;
            end
            a = [a l];
        end
    end
    a = sort(a);
end
