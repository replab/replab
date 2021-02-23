syms a b c d real;
s1 = [1 0; 0 1];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];
X = [(a+1i*b) (c+1i*d); (-c+1i*d) (a-1i*b)];

Y = kron(X, eye(2));

Z = [ a -b  c -d
      b  a  d  c
     -c -d  a  b
      d -c -b  a];

W = [s1-sy s1+sy; -s1-sy s1-sy]/2;
