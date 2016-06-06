function output = zlibdecode(input)
%ZLIBDECODE Decompress input bytes using ZLIB.
%
%    output = zlibdecode(input)
%
% The function takes a compressed byte array INPUT and returns inflated
% bytes OUTPUT. The INPUT is a result of GZIPENCODE function. The OUTPUT
% is always an 1-by-N uint8 array. JAVA must be enabled to use the function.
%
% See also zlibencode typecast

buffer = java.io.ByteArrayOutputStream();
zlib = java.util.zip.InflaterOutputStream(buffer);
zlib.write(input, 0, numel(input));
zlib.close();
output = typecast(buffer.toByteArray(), 'uint8')';

end
