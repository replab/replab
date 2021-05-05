function dts = testsInElement(el)
% Extracts all tests in a source element
%
% Args:
%   el (`+replab.+infra.SourceElement`): Source element to extract tests from
%
% Returns:
%   row cell array of `.DocTest`: Extracted doctests
    if el.doc.isempty
        dts = {};
        return
    end
    errFun = @(ln) replab.infra.doctests.errFunElement(el, el.doc.lineNumbers(ln));
    dts = replab.infra.doctests.parseTests(el.doc.lines, errFun);
    dts = cellfun(@(dt) dt.mapLineNumbers(@(i) el.doc.lineNumbers(i)), dts, 'uniform', 0);
end
