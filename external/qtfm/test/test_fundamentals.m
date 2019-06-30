function test_fundamentals
% Test code for the fundamental quaternion functions.

% Copyright ? 2005, 2006 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing quaternion fundamentals ...');

% check(qi == q1, 'qi and q1 are not equal.'); q1/q2/q3 are now obsolete.
% check(qj == q2, 'qj and q2 are not equal.');
% check(qk == q3, 'qk and q3 are not equal.');

check(qi .* qj .* qk == -1, 'ijk is not -1.');

check(i.^2 == -1, 'i^2 is not -1, perhaps been redefined?');

check(abs(qi) == 1, 'abs(qi) is not 1');
check(abs(qj) == 1, 'abs(qj) is not 1');
check(abs(qk) == 1, 'abs(qk) is not 1');

check(quaternion(1,1,1,1) == 1 + qi + qj + qk, 'Constructor error.');
check(abs(quaternion(1,1,1,1)) == 2, 'abs(1,1,1,1) is not 2.');
check(abs(quaternion(1,1,1,1)+quaternion(1,1,1,1).*i) == 2 + 2.*i,...
                                     'abs error on complex value.');

check(quaternion(5, 6, 7, 8) .* (1 + 2 .* qi + 3 .* qj + 4 .* qk) ...
    - quaternion(-60,20,14,32) == 0, 'Error in simple multiplication.');

check(quaternion(2,3,4,5) .* quaternion(5,6,7) - ...
      2 .* quaternion(5,6,7) - ...
      quaternion(3,4,5) .* quaternion(5,6,7) == 0, ...
      'Error in multiplication of full quaternion by pure.');
  
check(quaternion(5,6,7) .* quaternion(2,3,4,5) - ...
      2 .* quaternion(5,6,7) - ...
      quaternion(5,6,7) .* quaternion(3,4,5) == 0, ...
      'Error in multiplication of pure quaternion by full.');

check(quaternion(1,2,3).*2 - quaternion(2,4,6) == 0, ...
    'Error in multiplication of pure quaternion by numeric.');
check(2.*quaternion(1,2,3) - quaternion(2,4,6) == 0, ...
    'Error in multiplication of numeric by pure quaternion.');

check(quaternion(1,2,3,4).*2 - quaternion(2,4,6,8) == 0, ...
    'Error in multiplication of quaternion by numeric.');
check(2.*quaternion(1,2,3,4) - quaternion(2,4,6,8) == 0, ...
    'Error in multiplication of numeric by quaternion.');

check(                 quaternion(1,2,3).*quaternion(4,5,6) ...
      + scalar_product(quaternion(1,2,3), quaternion(4,5,6))...
      - vector_product(quaternion(1,2,3), quaternion(4,5,6)) == 0,...
      'Error in scalar/vector product test.');

check(scalar_product(quaternion(1,2,3,4), quaternion(5,6,7,8)) == 70,...
    'Error in scalar product test.');

compare(axis(quaternion(42,1,1,1)), axis(quaternion(3,3,3)), 1e-14,...
    'Error in axis comparison.');
                                                    
check(0 + [qi qj qk] == [0+qi 0+qj 0+qk], 'Addition fails test 1.');
check([qi qj qk] + 0 == [qi+0 qj+0 qk+0], 'Addition fails test 2.');
check(qi + (1 + qj)  == 1 + qi + qj,      'Addition fails test 3.');

check(quaternion(1,2,3,4) + quaternion(1,2,3)...
                         == quaternion(1,3,5,7),'Addition fails test 4.');
    
check(quaternion(1,2,3)   + quaternion(1,2,3,4)...
                         == quaternion(1,3,5,7),'Addition fails test 5.');
                     
check(quaternion(1,2,3,4) + quaternion(1,2,3,4)...
                         == quaternion(2,4,6,8),'Addition fails test 6.');
                     
check(quaternion(1,2,3)   + quaternion(1,2,3)...
                         == quaternion(2,4,6),'Addition fails test 7.');

check(conj(quaternion(1, 2, 3, 4)) == ...
           quaternion(1,-2,-3,-4), 'Conjugate fails test 1.');
check(conj(quaternion(1, 2, 3, 4), 'hamilton') == ...
           quaternion(1,-2,-3,-4), 'Conjugate fails test 2.');

check(conj(complex(quaternion(1,2,3,4), quaternion(5,6,7,8)), 'complex')...
       ==  complex(quaternion(1,2,3,4),-quaternion(5,6,7,8)), ...
                                   'Conjugate fails test 3.');
check(conj(complex(quaternion(1,2,3,4), quaternion(5,6,7,8)), 'total') == ...
           complex(quaternion(1,-2,-3,-4),-quaternion(5,-6,-7,-8)), ...
                                   'Conjugate fails test 4.');

disp('Passed');

% $Id: test_fundamentals.m 1004 2017-11-15 17:14:09Z sangwine $

