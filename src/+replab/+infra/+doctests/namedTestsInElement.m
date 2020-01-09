function s = namedTestsInElement(sel)
% Extracts all tests in a source element and its potential members
%
% When called on a class, it includes the tests of the class members.
%
% Args:
%   sel (`+replab.+infra.SourceElement`): Source element to extract tests from
%
% Returns:
%   struct: A structure whose field names are the test names, and values
%           are row cell vectors of `.DocTest`, with their line numbers
%           referring to the position in the original filename.
    s = struct;
    dt = replab.infra.doctests.testsInElement(sel);
    if ~isempty(dt)
        s.(sel.name) = dt;
    end
    if isa(sel, 'replab.infra.Class')
        cels = struct2cell(sel.ownElementsStruct);
        for i = 1:length(cels)
            cel = cels{i};
            dt = replab.infra.doctests.testsInElement(cel);
            if ~isempty(dt)
                if isfield(s, cel.name)
                    newValue = horzcat(s.(cel.name), dt);
                else
                    newValue = dt;
                end
                s.(cel.name) = newValue;
            end
        end
    end
end
