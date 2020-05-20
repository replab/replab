function elements = dimino(generators, addition, identity)

%%% INPUTS
%%% generators: a cell array of the generating elements
%%% addition: a function for adding two group elements
%%% identity: the identity of the group to be generated
%%%
%%% RETURNS elements: a cell array of all group elements 
%%%
%%% Ex. dimino({[2,3,1],[2,1,3]}, @S3.compose, [1,2,3]) returns all
%%% elements in the S3 permutation group

%%% Check generators are provided
if isempty(length(generators))
    return
end


%%% check if any generators are floats and find the max error
%%% check for differences in floats below error of the max value squared * 10
float = 0;
max_val = 0;
for i = 1:length(generators)
    if isfloat(generators{i})
        float = 1;
        max_gen = max(abs(generators{i}));
        max_gen = max_gen(1);
        if max_gen > max_val
            max_val = max_gen;
        end
    end
end
err = eps(max_val^2) * 10;

%%% Check that none of the generators are the identity
i = 1;
while i <= length(generators)
    if generators{i} == identity
        generators(i) = [];
    elseif float & abs(generators{i} - identity) < err
        generators(i) = [];
    end
    i = i+1;
end

% add check for whether generators are all from the same group?

%%% simple Dimino's algorithm

% add cyclic group of first generator
order = 1;
elements = {identity}; 
s1 = generators{1};
g = s1;
if float
    while abs(g - identity) > err
    order = order + 1;
    elements{order} = g;
    g = addition(g, s1);
    end
else
    while g ~= identity
        order = order + 1;
        elements{order} = g;
        g = addition(g, s1);
    end
end

for i = 2:length(generators)
    s = generators{i};
    not_redundant = 1;
    for e = 1:length(elements)
        if elements{e} == s
            not_redundant = 0;
            break
        elseif float & abs(elements{e} - s) < err % change this order somehow
            not_redundant = 0;
            break
        end
    end
    if not_redundant
        prev_order = order;
        order = order + 1;
        elements{order} = s;
        for j = 2:prev_order
            order = order + 1;
            elements{order} = addition(elements{j}, s);
        end

        cos_pos = prev_order + 1;
        while cos_pos <= order
            for n = 1:length(generators)
                sn = generators{n};
                elmt = addition(elements{cos_pos}, sn);
                add_coset = 1;
                for e = 1:length(elements)
                    if elements{e} == elmt
                        add_coset = 0;
                        break
                    elseif float & abs(elements{e} - elmt) < err
                        add_coset = 0;
                        break
                    end
                end
                if add_coset
                    order = order + 1;
                    elements{order} = elmt;
                    for j = 2:prev_order
                        order = order + 1;
                        elements{order} = addition(elements{j}, elmt);
                    end
                end
            end
            cos_pos = cos_pos + prev_order;
        end
        
    end
end

end