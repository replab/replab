function res = md5(str)
% Returns the MD5 hash of an ASCII string
%
% Args:
%   str (charstring): String to hash
%
% Returns:
%   charstring: MD5 hash in hexadecimal notation
    md = javaMethod('getInstance', 'java.security.MessageDigest', 'MD5');
    bytes = javaMethod('getBytes', javaObject('java.lang.String', str))
    res = typecast(md.digest(bytes), 'uint8');
    res = strjoin(arrayfun(@(b) dec2hex(b, 2), res, 'uniform', 0), '');
end
