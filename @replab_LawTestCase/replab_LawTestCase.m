function obj = replab_LawTestCase(lawsInstance, testName, lawMethodName, filePath, domains)
% Creates a test case for a law
%
% Args:
%   lawsInstance (`+replab.Laws`): Instance of the laws class
%   lawName (charstring): Name of the test
%   lawMethodName (charstring): Method in the laws class to test
%   filePath (charstring): Path to the laws class
%   domains (cell(1,\*) of `+replab.Domain`): List of domain elements to provide as arguments
    s = struct('lawsInstance', {lawsInstance}, 'lawMethodName', {lawMethodName}, 'domains', {domains});
    obj = class(s, 'replab_LawTestCase', MOxUnitTestCase(testName, filePath));
end
