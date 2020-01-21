function className = commonClass(groups)
% Returns the lowest denominator class type for the factor groups
    priority = 3;
    for i = 1:length(groups)
        group = groups{i};
        if isa(group, 'replab.NiceFiniteGroup')
            priority = min(priority, 3);
        elseif isa(group, 'replab.FiniteGroup')
            priority = min(priority, 2);
        elseif isa(group, 'replab.CompactGroup')
            priority = min(priority, 1);
        else
            error('Error, one of the elements is not a group');
        end
    end
    names = {'replab.CompactGroup' 'replab.FiniteGroup' 'replab.NiceFiniteGroup'};
    className = names{priority};
end
