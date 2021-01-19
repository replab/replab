function [xf, SS, cnt, res, XY] = LMFnlsq2(varargin)
% LMFNLSQ2   Solve one or more of (over)determined nonlinear equations
% in the least squares sense, prints Jacobian matrix and elapsed time.
%
% Obtained from https://www.mathworks.com/matlabcentral/fileexchange/39564-lmfnlsq2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A solution is obtained by a Fletcher's version of the Levenberg-Maquardt
% algoritm for minimization of a sum of squares of equation residuals.
% The main domain of LMFnlsq2 applications is in curve fitting during
% processing of experimental data (see Examle 3).
%
% [Xf, Ssq, CNT, Res, XY] = LMFnlsq2(FUN,Xo,Options)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Input arguments:
% FUN     is a function handle or a function M-file name (string) that
%         evaluates m-vector of equation residuals. Residuals are defined
%         as   res = FUN(x) - y,  where y is m-vector of given values.
% Xo      is n-vector of initial guesses of solution, unknown parameters,
%         n >= 1;
% Options is an optional set of Name/Value pairs of control parameters
%         of the algorithm. It may also be preset by calling:
%         Options = LMFnlsq2; %   for setting default values for Options,
%         Options = LMFnlsq2('default');  or by a set of Name/Value pairs:
%         Options = LMFnlsq2('Name',Value, ... ), or updating the Options
%                   set by calling
%         Options = LMFnlsq2(Options,'Name',Value, ...).
%
%    Name   Values {default}         Description
% 'Display'    {[0,0]}      Display control DC
%                             DC(1)= 0  no display
%                                   -k  don't display lambdas
%                                   |k| disp. initial & each k-th iteration
%                             DC(2)= 0  no display
%                                    1  disp. elapsed time one column of J
%                                    2  disp. columns of J matrix
% 'Printf'    name|handle   Function for displaying of results; {@printit}
% 'Jacobian'  name|handle   Jacobian matrix function; {@finjac}
% 'FunTol'      {1e-7}      norm(FUN(x),1) stopping tolerance;
% 'XTol'                    Stopping tolerance:
%               {1e-7}        scalar for equal delta x_i, or
%                             vector for particular delta x_i;
% 'Basdx'      {25e-9}      Basic step for Jacobian matrix evaluation
% 'MaxIter'     {100}       Maximum number of iterations;
% 'ScaleD'                  Scale control:
%               value         D = eye(m)*value;
%               vector        D = diag(vector);
%                {[]}         D(k,k) = JJ(k,k) for JJ(k,k)>0, or
%                                    = 1 otherwise,
%                                      where JJ = J.'*J
% 'Trace'        {0}        Tracing control:
%                             0 = don't save iteration results x^(k)
%                             nonzero = save iteration results x^(k)
% 'Lambda'       {0}        Initial value of parameter lambda
%
% A user may supply his own functions for building Jacobian matrix and for
% displaying intermediate results (corresponding to 'Jacobian' and
% 'Printf' options respectively).
% Not defined fields of the Options structure are filled by default values.
%
% Output Arguments:
%   Xf    final solution approximation
%   Ssq   sum of squares of residuals
%   Cnt     >0          count of iterations
%           -MaxIter,   no converge in MaxIter iterations
%   Res   number of calls of FUN and Jacobian matrix
%   XY    points of intermediate results in iterations
%
% Forms of function callings:
%   LMFnlsq2                             %   Display help to LMFnlsq2
%   Options = LMFnlsq2;                  %   Settings of Options
%   Options = LMFnlsq2('default');       %   The same as  Options = LMFnlsq2;
%   Options = LMFnlsq2(Name1,Value1,Name2,Value2,…);
%   Options = LMFnlsq2(Options,Name1,Value1,Name2,Value2,…);
%   x = LMFnlsq2(Eqns,x0);               %   Solution with default Options
%   x = LMFnlsq2(Eqns,x0,Options);       %   Solution with preset Options
%   x = LMFnlsq2(Eqns,x0,Name1,Value1,Name2,Value2,…);% W/O preset Options
%   [x,ssq] = LMFnlsq2(Eqns,x0,…);       %   with output of sum of squares
%   [x,ssq,cnt] = LMFnlsq2(Eqns,x0,…);   %   with iterations count
%   [x,ssq,cnt,nfJ,xy] = LMFnlsq2(Eqns,x0,…);
%
% Note:
% It is recommended to modify slightly the call of LMFnlsq2 for improvement
% of reliability of results and stability of the iteration process, namely
% if range of values of x's be large:
%   x0 = vector of user-defined initial guess;
%   x0 = x0(:);
%   x  = LMFnlsq2(Eqns,size(length(x0)), ...);
% This form has been used for testing problems of NIST. solution of one of
% them, BoxBOD, is included in the file BoxBOD.zip.
%
%% Example 1:
% The general Rosenbrock's function has the form of a sum of squares:
%    f(x) = 100(x(2)-x(1)^2)^2 + (1-x(1))^2
% Optimum solution gives f(x)=0 for x(1)=x(2)=1. Function f(x) can be
% expressed in the form
%    f(x) = f1(x)^2 =f2(x)^2,
% where f1(x) = 10(x(2)-x(1)^2) and f2(x) = 1-x(1).
% Values of the functions f1(x) and f2(x) can be used as residuals.
% LMFnlsq2 finds the solution of this problem in 17 iterations with default
% settings. The function FUN has the form of named function
%
%   function r = rosen(x)
%% ROSEN   Rosenbrock valey residuals
%   r = [ 10*(x(2)-x(1)^2)      %   first part,  f1(x)
%         1-x(1)                %   second part, f2(x)
%       ];
% or an anonymous function
%   rosen = @(x) [10*(x(2)-x(1)^2); 1-x(1)];
%
%% Example 2:
% Even that the function LMFnlsq2 is devoted for solving unconstrained
% problems, sometimes it is possible to solve also constrained problems,
% eg.:
% Find the least squares solution of the Rosenbrock's valey inside a circle
% of the unit diameter centered at the origin. In this case, it is necessary
% to build third function, a penalty, which is zero inside the circle and
% increasing outside it. This property has, say, the next function:
%    f3(x) = sqrt(x(1)^2 + x(2)^2) - d, where d is a distance from the
% circle border. Its implementation using named function has the form
%   function r = rosen(x)
%% ROSEN   Rosenbrock valey with a constraint
%   d = sqrt(x(1)^2+x(2)^2)-.5; %   distance from r=0.5
%   r = [ 10*(x(2)-x(1)^2)      %   first part,  f1(x)
%         1-x(1)                %   second part, f2(x)
%         (d>0)*d*w             %   penalty outside the feasible domain
%       ];                      %   w = 1000 is a weight of the condition
% or when anonymous functions are prefered
%    d     = @(x) sqrt(x'*x)-.5; %   A distance from the radius r=0.5
%    rosen = @(x) [10*(x(2)-x(1)^2); 1-x(1); (d(x)>0)*d(x)*1000];
%
% The solution is obtained in all cases by setting (say) FUN='rosen' in the
% first case, or FUN=rosen in the second one, and by calling
%    [x,ssq,cnt,loops]=LMFnlsq2(FUN,[-1.2,1],'Display',1,'MaxIter',50)
% yielding results
%    x=[0.45565; 0.20587],  |x|=0.5000,  ssq=0.29662,  cnt=43,  loops=58
% Exactly the same result is obtained, if analytical formula for the
% Jacobian matrix J is supplied. In case of anonymous function, the form of
% J may be defined as
%    jacob = @(x) [-20*x(1),10;-1,0;(d(x)>0)*1000/(d(x)+.5)*[x(1),x(2)]];
% and a minimum call is
%    [x,ssq,cnt,loop]=LMFnlsq2('rosen',[-1.2,1],'Jacobian',jacob);
% Since the formula for J is more complicated than that for the function in
% this case, the total time is shorter for evaluating J from finite
% differences, what is default.
%
%% Example 3 - Curve fit:
% For curve fit, see the accompanying script LMFnlsq2test and a document
% LMFnlsq2.pdf.
%
% Reference:
% Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least
% Squares. Rpt. AERE-R 6799, Harwell

% M. Balda,
% Institute of Thermomechanics,
% Academy of Sciences of The Czech Republic,
% balda AT cdm DOT cas DOT cz
% miroslav AT balda DOT cz
%
% 2007-07-02    a bit modified function LMFsolve
% 2007-10-08    Formal changes, improved description
% 2007-11-01    Completely reconstructed into LMFnlsq, new optional
%               parameters, improved stability
% 2007-12-06    Widened Option setting, improved help and description
% 2008-07-08    Complemented part for evaluation of Jacobian matrix from
%               an analytical formula. Small changes have been made in
%               description and comments.
% 2009-01-20    Implemented new subfunction for printing intermediate
%               results, and introduced a slight modification both of the
%               function and a testing script code. The pdf-file of
%               description has been complemented by a short theoretical
%               description of the LMF method.
%
% 2011-12-15 v.2.0 Improved description (help), improved printout
% 2012-02-27 v.2.1 Introduced new option 'Basdx' - vector of basic steps in
%                  iterations of x-solutions independent on XTol option.
%				   Renamed from LMFnlsq on LMFnlsq2
% 2012-12-05 v.2.2 Improved stability by introducing min(lambda)=1e-5
% 2013-04-02 v 2.3 Modified for the application of function LMFnlsq2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0 && nargout==0, help LMFnlsq2, return, end     %   Display help

%   Options = LMFnlsq2;
%   Options = LMFnlsq2('default');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Default Options
if nargin==0 || (nargin==1 && strcmpi('default',varargin(1)))
   xf.Display  = [0,0];     %   no print of iterations
   xf.Jacobian = 'finjac';  %   finite difference Jacobian approximation
   xf.MaxIter  = 0;         %   maximum number of iterations allowed
   xf.ScaleD   = [];        %   automatic scaling by D = diag(diag(J'*J))
   xf.FunTol   = 1e-7;      %   tolerace for final function value
   xf.XTol     = 1e-7;      %   tolerance on difference of x-solutions
   xf.Printf   = 'printit'; %   disply intermediate results
   xf.Trace    = 0;         %   don't save  intermediate results
   xf.Lambda   = 0;         %   start with Newton iteration
   xf.Basdx    = 25e-9;     %   basic step dx for gradient evaluation
   return

%   Options = LMFnlsq2(Options,name,value,...);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Updating Options
elseif isstruct(varargin{1}) % Options=LMFnlsq2(Options,'Name','Value',...)
    if ~isfield(varargin{1},'Jacobian')
        error('Options Structure not Correct for LMFnlsq2.')
    end
    xf=varargin{1};          %   Options
    for i=2:2:nargin-1
        name=varargin{i};    %   option to be updated
        if ~ischar(name)
            error('Parameter Names Must be Strings.')
        end
        name=lower(name(isletter(name)));
        value=varargin{i+1}; %   value of the option
        if strncmp(name,'d',1), xf.Display  = value;
        elseif strncmp(name,'f',1), xf.FunTol   = value;
        elseif strncmp(name,'x',1), xf.XTol     = value;
        elseif strncmp(name,'j',1), xf.Jacobian = value;
        elseif strncmp(name,'m',1), xf.MaxIter  = value;
        elseif strncmp(name,'s',1), xf.ScaleD   = value;
        elseif strncmp(name,'p',1), xf.Printf   = value;
        elseif strncmp(name,'t',1), xf.Trace    = value;
        elseif strncmp(name,'l',1), xf.Lambda   = value;
        elseif strncmp(name,'b',1), xf.Basdx    = value;
        else   disp(['Unknown Parameter Name --> ' name])
        end
    end
    return

%   Options = LMFnlsq2(name,value,...);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Pairs of Options
elseif ischar(varargin{1})  % check for Options=LMFnlsq2('Name',Value,...)
   Pnames=char('display','funtol','xtol','jacobian','maxiter','scaled',...
               'printf','trace','lambda');
   if strncmpi(varargin{1},Pnames,length(varargin{1}))
      xf=replab.numerical.LMFnlsq2('default');  % get default values
      xf=replab.numerical.LMFnlsq2(xf,varargin{:});
      return
   end
end

%   [Xf,Ssq,CNT,Res,XY] = LMFnlsq2(FUN,Xo,Options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               OPTIONS
%               *******
FUN=varargin{1};            %   function handle
if ~(isvarname(FUN) || isa(FUN,'function_handle'))
   error('FUN Must be a Function Handle or M-file Name.')
end
xc = varargin{2};           %   Xo
xc = xc(:);                 %   Xo is a column vector
if ~exist('options','var')
    options = replab.numerical.LMFnlsq2('default');
end
if nargin>2                 %   OPTIONS
    if isstruct(varargin{3})
        options=varargin{3};
    else
        for i=3:2:size(varargin,2)-1
            options=replab.numerical.LMFnlsq2(options, varargin{i},varargin{i+1});
        end
    end
else
    if ~exist('options','var')
        options = replab.numerical.LMFnlsq2('default');
    end
end




%               INITIATION OF SOLUTION
%               **********************
tic;
x = xc(:);
n = length(x);

bdx = options.Basdx;        %   basic step dx for gradient evaluation
lb  = length(bdx);
if lb==1
    bdx = bdx*ones(n,1);
elseif lb~=n
    error(['Dimensions of vector dx ',num2str(lb),'~=',num2str(n)]);
end

epsf  = options.FunTol(:);
ipr   = options.Display;
JAC   = options.Jacobian;
maxit = options.MaxIter;    %   maximum permitted number of iterations
if maxit==0, maxit=100*n; end
printf= options.Printf;

l  = options.Lambda;
lc = 1;
r  = feval(FUN,x);          %   initial "residuals"
%~~~~~~~~~~~~~~~~
SS = r'*r;

feval(printf,ipr,-1);       %   Table header
dx = zeros(n,1);
res= 1;
cnt=0;
feval(printf,ipr,cnt,res,SS,x,dx,l,lc) %    Initial state

[A,v] = getAv(FUN,JAC,x,r,bdx,ipr);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trcXY = options.Trace;      %   iteration tracing
if trcXY
    XY = zeros(n,maxit);
    XY(:,1) = x;
else
    XY = [];
end

D = options.ScaleD(:);      %   CONSTANT SCALE CONTROL D
if isempty(D)
    D = diag(A);            %   automatic scaling
else
    ld = length(D);
    if ld==1
        D = abs(D)*ones(n,1); %   scalar of unique scaling
    elseif ld~=n
        error(['wrong number of scales D, lD = ',num2str(ld)])
    end
end
D(D<=0)=1;
T = sqrt(D);

Rlo=0.25;
Rhi=0.75;

%               SOLUTION
%               ********    MAIN ITERATION CYCLE
while 1 %                   ********************
    if cnt>0
        feval(printf,ipr,cnt,res,SS,x,dx,l,lc)
    end
    cnt = cnt+1;
    if trcXY, XY(:,cnt+1)=x; end
    d = diag(A);
    s = zeros(n,1);
%                           INTERNAL CYCLE
    while 1 %               ~~~~~~~~~~~~~~
        while 1
            UA = triu(A,1);
%            l = max(l, 0.00001);    %   inserted on 2012-12-01
            A = UA'+UA+diag(d+l*D);
            [U,p] = chol(A);        %   Choleski decomposition
            %~~~~~~~~~~~~~~~
            if p==0, break, end
            l = 2*l;
            if l==0, l=1; end
        end
        dx = U\(U'\v);              %   vector of x increments
        vw = dx'*v;
        fin = -1;
        if vw<=0, break, end        %   The END

        for i=1:n
            z = d(i)*dx(i);
            if i>1, z=A(i,1:i-1)*dx(1:i-1)+z; end
            if i<n, z=A(i+1:n,i)'*dx(i+1:n)+z; end
            s(i) = 2*v(i)-z;
        end
        dq = s'*dx;
        s  = x-dx;
        rd = feval(FUN,s);
%            ~~~~~~~~~~~~
        res = res+1;
        SSP = rd'*rd;
        dS  = SS-SSP;
        fin = 1;
        if all((abs(dx)-bdx)<=0) || res>=maxit || abs(dS)<=epsf
            break                   %   The END
        end
        fin=0;
        if dS>=Rlo*dq, break, end
        A = U;
        y = .5;
        z = 2*vw-dS;
        if z>0, y=vw/z; end
        if y>.5, y=.5; end
        if y<.1, y=.1; end
        if l==0
            y = 2*y;
            for i = 1:n
                A(i,i) = 1/A(i,i);
            end
            for i = 2:n
                ii = i-1;
                for j= 1:ii
                    A(j,i) = -A(j,j:ii)*A(j:ii,i).*A(i,i);
                end
            end
            for i = 1:n
                for j= i:n
                    A(i,j) = abs(A(i,j:n)*A(j,j:n)');
                end
            end
            l  = 0;
            tr = diag(A)'*D;
            for i = 1:n
                z = A(1:i,i)'*T(1:i)+z;
                if i<n
                    ii = i+1;
                    z  = A(i,ii:n)*T(ii:n)+z;
                end
                z = z*T(i);
                if z>l, l=z; end
            end
            if tr<l, l=tr; end
            l  = 1/l;
            lc = l;
        end
        l = l/y;
        if dS>0, dS=-1e300; break, end
    end %  while            INTERNAL CYCLE LOOP
%                           ~~~~~~~~~~~~~~~~~~~

    if fin, break, end
    if dS>Rhi*dq
        l=l/2;
        if l<lc, l=0; end
    end
    SS=SSP;  x=s;  r=rd;
    [A,v] = getAv(FUN,JAC,x,r,bdx,ipr);
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % while                 MAIN ITERATION CYCLE LOOP
%                           *************************

if fin>0
    if dS>0
        SS = SSP;
        x  = s;
    end
end
if ipr(1)~=0
    disp(' ');
    feval(printf,sign(ipr),cnt,res,SS,x,dx,l,lc)
end
xf = x;
if trcXY, XY(:,cnt+2)=x; end
XY(:,cnt+3:end) = [];
if res>=maxit, cnt=-maxit; end
return
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,v] = getAv(FUN,JAC,x,r,bdx,ipr)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Calculate A, v, r
if isa(JAC,'function_handle')
    J = JAC(x);
else
    J = feval(JAC,FUN,r,x,bdx,ipr);
end
A = J'*J;
v = J'*r;
%end % getAv
% --------------------------------------------------------------------

function J = finjac(FUN,r,x,bdx,ipr)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  numer. approximation to Jacobi matrix
f =' %7.4f %7.4f %7.4f %7.4f';
rc = r(:);
lx = length(x);
J  = zeros(length(r),lx);
for k = 1:lx
    dx = bdx(k);
    xd = x;
    xd(k) = xd(k)+dx;
    rd = feval(FUN,xd);
%   ~~~~~~~~~~~~~~~~~~~
    J(:,k)=((rd(:)-rc)/dx);

% output columns of Jacobian matrix
  if ipr(1)~=0 && ipr(2)>0
    fprintf(' Column #%3d\n',k);
    tim  = toc;
    mins = floor(tim/60);
    secs = tim-mins*60;
    fprintf(' Elapsed time  =%4d min%5.1f sec\n',mins,secs)
    if ipr(2)==2
      fprintf([ f f f f '\n'], J(:,k)');
    end
    fprintf('\n');
  end
        %keyboard
        %interpt(k,dx,J,xd,rd);
end
%end % finjac
% --------------------------------------------------------------------

function printit(ipr,cnt,res,SS,x,dx,l,lc)
%        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  Printing of intermediate results
% For length(ipr)==1:
%  ipr(1) <  0  do not print lambda columns
%         =  0  do not print at all
%         >  0  print every (ipr)th iteration
%  ipr(2) =  0  do not prit time neither J
% For length(ipr)>1 && ipr(1)~=0:
%               print also elapsed time between displays
%  cnt = -1  print out the header
%        >0  print out results
if ipr(1)~=0
   if cnt<0                 %   table header
      disp('')
      nch = 50+(ipr(1)>0)*25;
      disp(char('*'*ones(1,nch)))
      fprintf('  itr  nfJ   SUM(r^2)        x           dx');
      if ipr(1)>0
          fprintf('           l           lc');
      end
      fprintf('\n');
      disp(char('*'*ones(1,nch)))
      disp('')
   else                     %   iteration output
      if rem(cnt,ipr(1))==0
          if ipr(2)>0
              tim  = toc;
              mins = floor(tim/60);
              secs = tim-mins*60;
              fprintf('\n Elapsed time  =%4d min%5.1f sec\n',mins,secs)
          end
          f='%12.4e ';
          if ipr(1)>0
             fprintf(['%5d %5d ' f f f f f '\n'],...
                 cnt,res,SS, x(1),dx(1),l,lc);
          else
             fprintf(['%5d %5d ' f f f '\n'],...
                 cnt,res,SS, x(1),dx(1));
          end
          for k=2:length(x)
             fprintf([blanks(25) f f '\n'],x(k),dx(k));
          end
      end
   end
end
%end % printit
%end % LMFnlsq2