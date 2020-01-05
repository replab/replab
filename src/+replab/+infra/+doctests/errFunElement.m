function errFunElement(el, ln)
    filename = el.absoluteFilename;
    lines = strsplit(fileread(filename), '\n', 'CollapseDelimiters', false);
    disp(['Parse error in ' replab.infra.linkOpen('%s:%d', '%s:%d', el.absoluteFilename, ln)]);
    disp(' ');
    disp(replab.infra.formatCodeContext(el, ln, 7));
    disp(' ');
end
