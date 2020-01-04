function dt = testsInElement(el)
% Extracts all tests in a source element
%
% Args:
%   el (`replab.infra.SourceElement`): Source element to extract tests from
%
% Returns:
%   row cell array of `.DocTest`: Extracted doctests
    if el.doc.isempty
        dt = {};
        return
    end
    ps = replab.infra.doctests.ParseState.fromDocTestBlock(el.doc.lines);
    
end
