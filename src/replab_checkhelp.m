function ok = replab_checkhelp
% Checks the soundness of documentation comments
%
% Returns:
%   logical: True if no problems were detected
    ok = true;
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
            ok = false;
            le = lasterror;
            last = p.consoleLine.lineContent;
            p.consoleLine.update('');
            if isa(el, 'replab.infra.SourceElement')
                replab.infra.doctests.errFunElement(el, el.startLineNumber);
                fprintf('\n');
            else
                fprintf('Parse error in inherited element %s\n', el.fullIdentifier);
                fprintf('\n');
            end
            fprintf('Error identifier: %s\n', le.identifier);
            disp(le.message);
            fprintf('\n\n');
            p.consoleLine.update(last);
            % Uncomment below for stack display
            % for i = 1:length(le.stack)
            %     disp(strjoin(replab.longStr(le.stack(i)), '\n'));
            % end
        end
    end
    p.finish('Check done');
end
