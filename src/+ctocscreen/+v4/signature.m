function s=signature()
%SIGNATURE Include all project MATLAB dependencies, including private Lambert.
folder=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
files=dir(fullfile(folder,'src','**','*.m'));
extra=dir(fullfile(folder,'scripts','*v4*.m')); files=[files;extra];
paths=arrayfun(@(f)fullfile(f.folder,f.name),files,'UniformOutput',false);
[paths,order]=sort(paths); files=files(order); %#ok<NASGU>
s=struct('interface','free_trajectory_v4_8','paths',{{}},'sha256',{{}});
for k=1:numel(paths)
 s.paths{k}=strrep(paths{k},[folder filesep],''); s.sha256{k}=ctocscreen.v3FileHash(paths{k});
end
s.target_csv_sha256=ctocscreen.v3FileHash(fullfile(folder,'data','ctoc14b_targets.csv'));
end
