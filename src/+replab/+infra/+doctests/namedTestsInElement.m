function s = namedTestsInElement(el)
% Extracts all tests in a source element and its potential members
%
% When called on a class, it includes the tests of the class members.
%
% Args:
%   el (`replab.infra.SourceElement`): Source element to extract tests from
%
% Returns:
%   struct: A structure whose field names are the test names, and values
%           are row cell vectors of `.DocTest`, with their line numbers
%           referring to the position in the original filename.
    s = struct;
    doctests = replab.infra.DocTest.
end
