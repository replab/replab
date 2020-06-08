function classes = ConjugacyClasses(elements, generators, group)
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


classes = {};
i = 0;

for d = 1:length(elements)
    curr_elt = elements{d};
    % Check whether the element is already in a conjugacy class
    not_in_class = 1;
    for j = 1:length(classes)
        for k = 1:length(classes{j})
            if group.eqv(classes{j}{k}, curr_elt)
                not_in_class = 0;
                break
            end
        end
        if k < length(classes{j})
            break
        end
    end
    if not_in_class
        % create a new class array if element is not already in a class
        i = i + 1;
        new_class = elements(d);
        pos = 2; % next position to add element to class
        ci = 1; 
        % continue until conjugates of every element in set have been added
        while ci <= length(new_class)
            g = new_class{ci};
            % conjugate elements with respect to all generators
            for j = 1:length(generators)
                s = generators{j};
                conjugate = group.compose(group.inverse(s), ...
                    group.compose(g, s));
                add_to_class = 1;
                % check whether conjugate is in the class already
                for k = 1:length(new_class)
                    if group.eqv(conjugate, new_class{k})
                        add_to_class = 0;
                        break
                    end
                end
                if add_to_class
                    new_class{pos} = conjugate;
                    pos = pos + 1;
                end
            end
            ci = ci + 1; % conjugate element following g
        end
        classes{i} = new_class; % start next class
    end
end

end