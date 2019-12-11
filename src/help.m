function help(name)
    persistent codeBase
    if isequal(name(1:6), 'replab')
        if isempty(codeBase)
            [srcRoot, name, ~] = fileparts(mfilename('fullpath'));
            codeBase = replab.infra.CodeBase.crawl(srcRoot);
        end
        disp(codeBase.lookupName(name))
    else
        error('JDB: put your special code here');
    end
end
