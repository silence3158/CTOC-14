function hash=implementationHash()
%IMPLEMENTATIONHASH Bind checkpoints to the actual MATLAB package sources.
folder=fileparts(mfilename('fullpath')); files=dir(fullfile(folder,'*.m'));
[~,idx]=sort({files.name}); files=files(idx);
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(files)
 md.update(typecast(unicode2native(files(k).name,'UTF-8'),'int8'));
 fid=fopen(fullfile(folder,files(k).name),'rb'); cleanup=onCleanup(@()fclose(fid));
 bytes=fread(fid,Inf,'*uint8'); md.update(typecast(bytes,'int8')); clear cleanup;
end
bytes=typecast(md.digest(),'uint8'); hash=lower(reshape(dec2hex(bytes,2).',1,[]));
end
