function h = v3FileHash(filename)
%V3FILEHASH SHA-256 of exact bytes, including non-ASCII paths.
fid=fopen(filename,'rb');
assert(fid>=0,'ctocscreen:v3:missingSource','Cannot open %s',filename);
c=onCleanup(@()fclose(fid));
b=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');
h=lower(reshape(dec2hex(typecast(md.digest(typecast(b,'int8')),'uint8'),2).',1,[]));
end
