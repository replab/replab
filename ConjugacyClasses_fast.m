function class_list = ConjugacyClasses_fast(elements, generators, group)
% generates the conjugacy classes all elements of a group
% 
% Orbit computation algorithm:
% Butler, Gregory. “Orbits and Schreier Vectors.” FUNDAMENTAL ALGORITHMS FOR PERMUTATION GROUPS, 
% by Gregory Butler, Springer-Verlag, 1991, pp. 57-59.
% 
% Args:
%   elements (cell array): elements of the group
%   generators (cell array): generating elements of the group
%   group (replab Group object): group used for composition, equivalence,
%                                   and identity
%
% Returns 
%   classes (cell array or cell arrays): all conjugacy classes 
%

%%% Determine group type
id = group.identity;
if iscell(id) % direct product
    vec = 0; mat = 0; dir_prod = 1;
    vecid = vectorize(id);
    elt_len = length(vecid);
elseif isvector(id) % permutation
    vec = 1; mat = 0; dir_prod = 0;
    elt_len = length(id);
elseif ismatrix(id) % integer matrix
    vec = 0; mat = 1; dir_prod = 0;
    dim = size(id);
    elt_len = numel(id);
else
    err('Error: data type not accepted')
end

%%% Store inverses of generators for conjugation
inverses = cell(1, length(generators));
for i = 1:length(generators)
    inverses{i} = group.inverse(generators{i});
end

classes = {};
i = 0;
approx_class_elmts = length(elements)/2; % probably can make this smaller

for d = 1:length(elements)
    if vec
        curr_elmt = elements{d};
    elseif mat
        curr_elmt = reshape(elements{d}, [1, elt_len]);
    else
        curr_elmt = vectorize(elements{d});
    end
    % Check whether the element is already in a conjugacy class
    not_in_class = 1;
    for j = 1:length(classes)
        if any(all(classes{j} - curr_elmt == 0, 2))
            not_in_class = 0;
            break
        end
    end
    if not_in_class
        % create a new class array if element is not already in a class
        i = i + 1;
        new_class = zeros(approx_class_elmts, elt_len);
        if vec
            new_class(1, :) = elements{d};
        elseif mat
            new_class(1, :) = reshape(elements{d}, [1, elt_len]);
        else
            new_class(1, :) = vectorize(elements{d});
            class_arr = elements(d);
        end
        pos = 1; % current position in class list
        ci = 1; 
        % continue until conjugates of every element in set have been added
        while ci <= pos
            if vec
                g = new_class(ci, :);
            elseif mat
                g = reshape(new_class(ci, :), dim);
            else
                g = class_arr{ci};
            end
            % conjugate elements with respect to all generators
            for j = 1:length(generators)
                s = generators{j};
                sinv = inverses{j};
                conjugate = group.compose(sinv, ...
                    group.compose(g, s));
                if vec
                    conj_vec = conjugate;
                elseif mat
                    conj_vec = reshape(conjugate, [1, elt_len]);
                else
                    conj_vec = vectorize(conjugate);
                end
                % add to class if not already present 
                if ~any(all(new_class - conj_vec == 0, 2))
                    new_class(pos + 1, :) = conj_vec;
                    if dir_prod
                        class_arr{pos + 1} = conjugate;
                    end
                    pos = pos + 1;
                end
            end
            ci = ci + 1; % conjugate element following g
        end
        % add new_class to class list
        classes{i} = new_class(1:pos, :);
        if dir_prod
            class_list{i} = class_arr;
        end
%         if vec
%             classes{i} = mat2cell(elements, ones(1, pos));
%         elseif mat
%             new_cell_arr = mat2cell(elements, ones(1, pos));
%             for i = 1:length(new_cell_arr)
%                 new_cell_arr{i} = reshape(new_cell_arr{i}, dim);
%             end
%             classes{i} = new_cell_arr;
%         else
%             classes{i} = class_arr; 
%         end
    end
end

if ~dir_prod
    class_list = cell(1, length(classes));
    for i = 1:length(classes)
        if vec
            mat_dim = size(classes{i});
            class_list{i} = mat2cell(classes{i}, ones(1, mat_dim(1)));
        else
            mat_dim = size(classes{i});
            new_cell_arr = mat2cell(classes{i}, ones(1, mat_dim(1)));
            for j = 1:length(new_cell_arr)
                new_cell_arr{j} = reshape(new_cell_arr{j}, dim);
            end
            class_list{i} = new_cell_arr;
        end
    end
end

function vector = vectorize(elmt)
    cell_locs = cellfun(@iscell, elmt);
    while any(cell_locs)
        elmt = [elmt(~cell_locs), elmt{cell_locs}];
        cell_locs = cellfun(@iscell, elmt);
    end
    vector = cell2mat(elmt);
end

end