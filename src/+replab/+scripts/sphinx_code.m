function sphinx_code(codeBase, targetFolder)
% Preprocesses source code to complement the Sphinx Matlab domain
%
% Preprocesses the source code of the given code base to complement the
% `Sphinx Matlab domain <https://github.com/sphinx-contrib/matlabdomain>`_
% job. Processed source files are written into the desired target folder.
% For all processed files, this function generates method and property
% occurences for inherited members and shifts the documentation of
% properties to be explicit (to support multiline documentation and types).
% This function also generates a table of contents for all source code
% files.
%
% This function erases any file in ``targetFolder``.
%
% Args:
%   codeBase (`+replab.infra.CodeBase`): Source code to be processed
%   targetFolder (charstring): Location for the processed code

    logFun = @(str) disp(str);

    % Clean the target folder
    [base, name] = fileparts(targetFolder);
    replab.infra.mkCleanDir(base, name, logFun);

    % Generate preprocessed source files for Sphinx
    logFun('Generating rich source code');
    els = codeBase.allSourceElements;
    pb = replab.infra.repl.ProgressBar(length(els));
    for i = 1:length(els)
        pb.step(i, els{i}.fullIdentifier);
        replab.infra.sphinx.writeEnrichedSource(targetFolder, els{i});
    end
    pb.finish;
    tocFile = fullfile(targetFolder, 'index.rst');
    replab.infra.sphinx.writeTOC(tocFile, codeBase);
end
