function result = matlabLaunchTestCovering
% This function creates a code covering report on matlab
%
% Note: this file is meant to be run from the root 'replab' folder

% This only works on matlab
assert(~replab.compat.isOctave);

import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.CodeCoveragePlugin;
import matlab.unittest.plugins.codecoverage.CoberturaFormat;


suite = TestSuite.fromFile('tests/matlabTestSuite.m');
runner = TestRunner.withTextOutput;
reportFile = 'coverage.xml';
reportFormat = CoberturaFormat(reportFile);

runner.addPlugin(CodeCoveragePlugin.forFolder([pwd, '/src'], 'IncludingSubfolders', true, 'Producing', reportFormat));
report = runner.run(suite);
result = report.Passed;
    
end
