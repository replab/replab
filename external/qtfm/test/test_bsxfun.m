function test_bsxfun
% Test code for the quaternion/octonion bsxfun functions.

% Copyright (c) 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing bsxfun ...')

% Method: compare quaternion bsxfun against Matlab bsxfun with the same
% data in the scalar part of the quaternion. We are not testing whether the
% function called by bsxfun works, but whether bsxfun correctly implements
% singleton expansion.

test_one_case([2,1], [1,2]); % Vector test cases.
test_one_case([3,1], [1,3]);
test_one_case([5,1], [1,4]);

test_one_case([2,2], [1,2]); % Matrix test cases.
test_one_case([3,4], [1,4]);
test_one_case([5,3], [1,3]);;

test_one_case([2,2,2], [1,2,2]);
test_one_case([2,2,2], [2,1,2]);
test_one_case([2,2,2], [2,2,1]);
test_one_case([10,10,1], [10,1,10]);

test_one_case([1,2,2,2], [2,2,2,2]);
test_one_case([2,1,2,2], [2,2,2,2]);
test_one_case([2,2,1,2], [2,2,2,2]);
test_one_case([2,2,2,1], [2,2,2,2]);
test_one_case([2,1,2,1], [2,2,2,2]);
test_one_case([1,2,1,1], [2,2,2,2]);

disp('Passed');

end

function test_one_case(v1, v2)
% v1 and v2 are two vectors of sizes. We use these to make arrays then
% apply bsxfun to the two arrays.

A = randi(100, v1); % Use integer test data, so we can use check and not
B = randi(100, v2); % worry about tolerances

C = bsxfun(@plus, A, B); % This is the Matlab bsxfun.
D = bsxfun(@plus, quaternion(A), quaternion(B)); % The QTFM bsxfun.
E = bsxfun(@plus,   octonion(A),   octonion(B)); % The QTFM bsxfun.

check(C == D, ['quaternions bsxfun fails test with sizes ', ...
                num2str(v1), ' and ', num2str(v2)]);
check(C == E, ['octonion bsxfun fails test with sizes ', ...
                num2str(v1), ' and ', num2str(v2)]);
            
A = randi(100, v2); % Now try with the sizes reversed.
B = randi(100, v1);

C = bsxfun(@times, A, B); % This is the Matlab bsxfun, this time we use .*
D = bsxfun(@times, quaternion(A), quaternion(B)); % The QTFM bsxfun.
E = bsxfun(@times,   octonion(A),   octonion(B)); % The QTFM bsxfun.

check(C == D, ['quaternion bsxfun fails test with sizes ', ...
                num2str(v2), ' and ', num2str(v1)]);
check(C == E, ['octonion bsxfun fails test with sizes ', ...
                num2str(v2), ' and ', num2str(v1)]);

end

% $Id: test_bsxfun.m 1004 2017-11-15 17:14:09Z sangwine $
