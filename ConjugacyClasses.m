function classes = ConjugacyClasses(elements, generators, group)
classes = {};
i = 0;
for d = 1:length(elements)
    curr_elt = elements{d};
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
        i = i + 1;
        new_class = elements(d);
        pos = 2;
        ci = 1;
        while ci <= length(new_class)
            g = new_class{ci};
            for j = 1:length(generators)
                s = generators{j};
                conjugate = group.compose(group.inverse(s), ...
                    group.compose(g, s));
                add_to_class = 1;
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
            ci = ci + 1;
        end
        classes{i} = new_class;
    end
end

end