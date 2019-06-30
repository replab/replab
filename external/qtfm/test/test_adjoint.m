function test_adjoint
% Test code for the quaternion adjoint and unadjoint functions.

% Copyright (c) 2007, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

disp('Testing adjoint/unadjoint ...');

% Create test matrices in various sizes.

for N=1:9

    q = randq(N);
    
    if any(any(unadjoint(adjoint(q)) ~= q))
        error(['quaternion/adjoint/unadjoint failed test 1 ',...
               'with matrix of size ', num2str(N)]);
    end
    
    if any(any(unadjoint(adjoint(q, 'real'), 'real') ~= q))
        error(['quaternion/adjoint/unadjoint failed test 2 ',...
               'with matrix of size ', num2str(N)]);
    end
    
    if any(any(unadjoint(adjoint(q, 'block'), 'block') ~= q))
        error(['quaternion/adjoint/unadjoint failed test 3 ',...
               'with matrix of size ', num2str(N)]);
    end
    
    if any(any(unadjoint(adjoint(q, 'real', 'block'), ...
                                    'real', 'block') ~= q))
        error(['quaternion/adjoint/unadjoint failed test 4 ',...
               'with matrix of size ', num2str(N)]);
    end

end

% Create a complex quaternion test matrix.

q = complex(randq(4), randq(4));

% Test 5 is not exact, like the other tests, so we must use compare and
% give a tolerance for the comparison. This is because the complex adjoint
% of a complex quaternion matrix mixes the elements of the original matrix
% in making the adjoint, whereas in other cases, the elements are simply
% repositioned with or without negation yielding exact reconstruction.

T = 1e-12;

compare(q, unadjoint(adjoint(q, 'complex'), 'complex'), T, ...
    'quaternion/adjoint/unadjoint failed test 5.');

if any(any(unadjoint(adjoint(q, 'quaternion'), 'quaternion') ~= q))
    error('quaternion/adjoint/unadjoint failed test 6.');
end

disp('Passed');

% $Id: test_adjoint.m 1004 2017-11-15 17:14:09Z sangwine $
