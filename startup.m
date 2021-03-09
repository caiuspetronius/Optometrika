function startup
%STARTUP. Adds the current directory to path.

pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(genpath(pwd), pathCell));
else
  onPath = any(strcmp(genpath(pwd), pathCell));
end

if (~onPath)
    addpath(genpath(pwd));
end

end

