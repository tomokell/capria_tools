% From https://uk.mathworks.com/matlabcentral/answers/56363-fully-resolving-path-names
% Downloaded 2/11/21
function absPath = getAbsPath(obj)
getAbsFolderPath = @(y) string(unique(arrayfun(@(x) x.folder, dir(y), 'UniformOutput', false)));
getAbsFilePath = @(y) string(arrayfun(@(x) fullfile(x.folder, x.name), dir(y), 'UniformOutput', false));
if isfolder(obj)
    absPath = getAbsFolderPath(obj);
elseif isfile(obj)
    absPath = getAbsFilePath(obj);
else
    error('The specified object does not exist.');
end