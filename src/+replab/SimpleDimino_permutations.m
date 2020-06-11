function elements = SimpleDimino_permutations(group)
% Returns all elements of the permutation group
% 
% Dimino's algorithm:
% Butler, Gregory. “List of Elements.” FUNDAMENTAL ALGORITHMS FOR PERMUTATION GROUPS, 
% by Gregory Butler, Springer-Verlag, 1991, pp. 14–23.
% 
% Args:
%   group (replab Permutation group): group used for composition, equivalence,
%                                     identity, and generators
%
% Returns 
%   elements (cell array): all group elements 
%
% Example:
% S3 = replab.S(3)
% replab.SimpleDimino({[2,3,1],[2,1,3]}, S3) returns all elements in the 
%                                               S3 permutation group


%%% Check generators are provided
generators = group.generators;

%%% simple Dimino's algorithm

% add cyclic group of first generator
order = 1;
elements = {group.identity};
tab = [Index(group.identity)];
s1 = generators{1};
g = s1;
% add powers of s1 until getting back to identity
while ~group.eqv(g, group.identity) 
    order = order + 1;
    elements{order} = g;
    tab = [tab, Index(g)];
    g = group.compose(g, s1);
end
tab = sort(tab);

% inductive steps in algorithm
for i = 2:length(generators)
    s = generators{i};
    index = Index(s);
    [redundant, tab] = Indexing(index, tab);
    if ~redundant
        prev_order = order;
        order = order + 1;
        elements{order} = s;
        % compose previous elements with new generator
        for j = 2:prev_order
            order = order + 1;
            new_elmt = group.compose(elements{j}, s);
            elements{order} = new_elmt;
            tab = sort([tab, Index(new_elmt)]);
        end

        % use coset to add new elements
        cos_pos = prev_order + 1; % first element in coset is generator s
        while cos_pos <= order
            for n = 1:length(generators)
                sn = generators{n};
                % left multiply the elements of the coset
                elmt = group.compose(elements{cos_pos}, sn);
                index = Index(elmt);
                [not_coset, tab] = Indexing(index, tab);
                if ~not_coset % add elmt to coset if not in elements
                    order = order + 1;
                    elements{order} = elmt;
                    tab = sort([tab, Index(elmt)]);
                    % right multiply all previous elements with the new
                    % element in the coset
                    for j = 2:prev_order
                        order = order + 1;
                        new_elmt = group.compose(elements{j}, elmt);
                        elements{order} = new_elmt;
                        tab = sort([tab, Index(new_elmt)]);
                    end
                end
            end
            % move the coset position to the position of elmt
            cos_pos = cos_pos + prev_order;
        end
        
    end
end

end
