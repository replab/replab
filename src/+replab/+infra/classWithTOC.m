function str = classWithTOC(codeBase, class)
% Returns the source code of a class with an injected Table of Contents
%
% Args:
%   class (`.Class`): Class to process
%
% Returns:
%   charstring: The processed source code
    lines = strsplit(fileread(class.fullFilename), '\n', 'CollapseDelimiters', false);
    pos = 2;
    while pos <= length(lines)
        l = lines{pos};
        if isempty(l) || l(1) ~= '%'
            break
        end
        pos = pos + 1;
    end
    t = cell(1, 3);
    t{1} = '% .. raw:: html';
    t{2} = '%    <h3>Methods</h3>';
    t{3} = '%';
    names = class.inheritedMemberNames(codeBase);
    for i = 1:length(names)
        name = names{i};
        if class.hasMember(name)
            member = class.member(name);
            t{end+1} = sprintf('%% - `%s` -- %s', name, member.doc.firstLine);
        else
            members = class.findInheritedMember(codeBase, name);
            member = members{1};
            t{end+1} = sprintf('%% - `%s` -- %s', member.sphinxFullName, member.doc.firstLine);
        end
    end
    t{end+1} = '%';
    t{end+1} = '% .. raw:: html';
    t{end+1} = '%';
    t{end+1} = '%    <hr>';
    str = strjoin(horzcat(lines(1:pos-1), t, lines(pos:end)), char(10));
end
