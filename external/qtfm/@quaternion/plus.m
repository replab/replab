function q = plus(l, r)
% +   Plus.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2008, 2010, 2016 Stephen J. Sangwine
%                                      and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Three cases have to be handled:
%
% l is a quaternion, r is not,
% r is a quaternion, l is not,
% l and r are both quaternions.

% An additional complication is that the parameters may be empty (numeric,
% one only) or empty (quaternion, one or both). Matlab adds empty + scalar
% to give empty: [] + 2 gives [], not 2, but raises an error (Matrix
% dimensions must agree) on attempting to add empty to an array. This
% behaviour is copied here, whether the empty parameter is a quaternion or
% a numeric empty. The result is always a quaternion empty.

ql = isa(l, 'quaternion');
qr = isa(r, 'quaternion');

% Find the sizes and classes of the arguments or their components if they
% are quaternions. We avoid calling the quaternion size function because of
% the overhead. Instead we can directly compute the size using the x
% component, which always exists.

if ql
    ls = size(l.x); lc = class(l.x);
else
    ls = size(l);   lc = class(l);
end

if qr
    rs = size(r.x); rc = class(r.x);
else
    rs = size(r);   rc = class(r);
end

if ~strcmp(lc, rc)
    % TODO Eliminate this error, if the non-matching class is logical. We
    % can handle this. At present qi + true fails with this error, whereas
    % qi - true will work because uminus converts true to double. This is a
    % stupid inconsistency.
    error(['Class of left and right parameters does not match: ', ...
        lc, ' ', rc])
end

pl = prod(ls); % The product of the elements of size(l).
pr = prod(rs); % The product of the elements of size(r).

sl = pl == 1; % Equivalent to isscalar(l), since ls must be [1, 1]. Note 1.
sr = pr == 1; % Equivalent to isscalar(r).

el = pl == 0; % Equivalent to isempty(l), since ls must be [0, 0], without
er = pr == 0; % calling an isempty function (quaternion or Matlab). Note 2.

% Now check for the case where one or other argument is empty and the other
% is scalar. In this case we choose to return an empty quaternion, similar
% to the behaviour of Matlab's plus function.

if (el && sr) || (sl && er)
    q = quaternion.empty;
    return
end

% Having now eliminated the cases where one parameter is empty and the
% other is scalar, the parameters must now either match in size, or one
% must be scalar. The concept of 'matching size' changed with Matlab R2016b
% as from that release Matlab accepts singleton dimensions as matching
% non-singleton dimensions (implicit singleton expansion). Previously there
% were checks here on the sizes, but the complexity of doing this with
% implicit expansion made it simpler to remove the checks and allow Matlab
% to perform them when the built-in plus function is called.
%
% It could be that both parameters are empty arrays of size [0, n] or
% [n, 0] etc, in which case we must return an empty array of the same size
% to match Matlab's behaviour.

if el && er
   % We must return an empty (quaternion) matrix.
   if ql && qr
       q = l; % r would do just as well, since both are quaternions.
   elseif ql
       q = l; % This must be l because r isn't a quaternion.
   else
       q = r; % This must be r because l isn't a quaternion.
   end
   return
end

if ql && qr

    % The scalar part(s) could be empty, and we have to handle this
    % specially, because [] + x is [] and not x, for some reason, in
    % Matlab, as noted above, whereas clearly if we add a pure quaternion
    % to a full quaternion, we want the result to have the scalar part of
    % the full quaternion and not be empty. In order to implement singleton
    % expansion we need to add an explicit array of zeros matching the size
    % of the vector parts so that Matlab's plus function does the singleton
    % expansion for us.
    
    q = l; % Copy one of the input parameters to avoid a constructor.
    
    q.x = q.x + r.x; % Deal with the vector part first, so we can return if
    q.y = q.y + r.y; % both scalar parts are empty.
    q.z = q.z + r.z;
    
    if isempty(l.w)
       if isempty(r.w)
           % Both scalar parts are empty. Since q.w is empty because we
           % copied it from l, we are done.
           return
       else
           % l.w is empty, but r.w is not.
           q.w = zeros(ls, lc) + r.w;
       end
    elseif isempty(r.w)
        % r.w is empty, but l.w is not.
        q.w = l.w + zeros(rs, rc);
    else % Neither l.w, nor r.w is empty.
        q.w = l.w + r.w;
    end


elseif isa(r, 'numeric')
        
    % The left parameter must be a quaternion, otherwise this function
    % would not have been called. We handle this by promoting the right
    % parameter to a quaternion, then adding it. This is not efficient, but
    % it is the only sensible way, for the moment, to provide implicit
    % singleton expansion. To see why, consider the case A + B where A is a
    % quaternion array of dimension 2 by 2, and B is a double of dimension
    % 2 by 1 by 2. Adding B to the scalar part of A is easy (Matlab will do
    % the singleton expansion), but how do we expand the vector part of A
    % to match? Converting B to quaternion and then adding solves this.

    q = l + quaternion(r); % TODO Eliminate the recursive call here?
    
elseif isa(l, 'numeric')
    
    q = quaternion(l) + r;

else
  error('Unhandled parameter types in function +/plus')
end

% Note 1. Actually ls could be [1,1,1] or [1,1,1,...], but these cases are
% treated as scalar by Matlab (try it with rand(1,1,1) to see).

% Note 2. Actually ls could be [0,1,2,....], which is also empty. However,
% Matlab will add two arrays of this type of the same size to yield another
% empty array of the same size.

% $Id: plus.m 1004 2017-11-15 17:14:09Z sangwine $
