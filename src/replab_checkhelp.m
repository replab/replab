function replab_checkhelp
% Checks the soundness of documentation comments
    help('--clear');

    rp = replab.settings.replabPath;
    srcRoot = fullfile(rp, 'src');


    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);

    disp('Checking');
    ac = cb.allClasses;
    elements = cb.allFunctions;
    for i = 1:length(ac)
        c = ac{i};
        elements = horzcat(elements, {c}, c.allElements);
    end

    p = replab.infra.repl.ProgressBar(length(elements));
    for i = 1:length(elements)
        el = elements{i};
        p.step(i, el.fullIdentifier);
        try
            evalc(sprintf('help %s', el.fullIdentifier));
            evalc(sprintf('help -f %s', el.fullIdentifier));
        catch
            le = lasterror;
            fprintf('\n');
            fprintf('Error in documentation of %s (line %d)\n', el.fullIdentifier, el.startLineNumber);
            fprintf('Error identifier: %s\n', le.identifier);
            disp(le.message);
            fprintf('\n\n');
            %for i = 1:length(le.stack)
            %    disp(strjoin(replab.longStr(le.stack(i)), '\n'));
            %end
        end
    end
    p.finish('Check done');
end
