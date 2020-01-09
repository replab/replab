function replab_checkhelp
% Checks the soundness of documentation comments
    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');

    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);

    disp('Checking functions');
    af = cb.allFunctions;
    for i = 1:length(af)
        f = af{i};
        fprintf('Checking %s\n', f.fullIdentifier);
        evalc(sprintf('help %s', f.fullIdentifier));
        evalc(sprintf('help -f %s', f.fullIdentifier));
    end
    ac = cb.allClasses;
    for i = 1:length(ac)
        c = ac{i};
        fprintf('Checking %s\n', f.fullIdentifier);
        evalc(sprintf('help %s', f.fullIdentifier));
        evalc(sprintf('help -f %s', f.fullIdentifier));
        els = c.allElements;
        for j = 1:length(els)
            el = els{j};
            fprintf('Checking %s\n', el.fullIdentifier);
            evalc(sprintf('help %s', el.fullIdentifier));
            evalc(sprintf('help -f %s', el.fullIdentifier));
        end
    end
end
