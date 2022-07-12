% Simple function to pass the input arguments stored as a struct in opts,
% using default values stored in the struct defaults where not specified
% within opts.
%
% Tom Okell, June 2022
%
% out = topassargs(opts,defaults,printdefaults,errorwithunknownopts,prefix)
%
% If errorwithunknownopts is true, an error will be thrown if a field of
% opts is not present in the defaults. prefix is used to be able to call
% this recursively when a struct contains a field which is another struct.
function out = topassargs(opts,defaults,printdefaults,errorwithunknownopts,prefix)

if nargin < 3; printdefaults = false; end
if nargin < 4; errorwithunknownopts = true; end
if nargin < 5; prefix = ''; end

% Initialise the output with the defaults
out = defaults;

% If passed an empty input, just return the defaults, which have already
% been assigned
if isempty(opts)
    disp('Using default values for all options')
    disp(out)
    return
end

% Assign the value in opts if present
params = fieldnames(defaults);
NoParamsAssigned = 0;
for ii = 1:length(params)
    
    % If the field is a struct, call the function recursively
    if isstruct(eval(['out.' params{ii}]))
        
        if ~isfield(opts, params{ii}) % If there is no field, pass an empty value
            eval(['opts.' params{ii} ' = [];'])
        end
        if isempty(prefix)
            nextprefix = [params{ii} '.'];
        else
            nextprefix = [prefix '.' params{ii} '.'];
        end
        eval(['out.' params{ii} ' = topassargs(opts.' params{ii} ', out.' params{ii} ',printdefaults,errorwithunknownopts,nextprefix);']);
    
    else % Not a struct option
        if isfield(opts, params{ii}) % If the option is present in the input, assign it
            eval(['out.' params{ii} ' = opts.' params{ii} ';'])
            NoParamsAssigned = NoParamsAssigned + 1;
        
        else % Need to pass in the default value
            
            defaultval = eval(['defaults.' params{ii}]);
            if (size(defaultval,1) > 1) || (length(defaultval)>50)
                str = 'Too big to display';
            else
                str = num2str(defaultval);
            end
            
            if printdefaults
                disp(['Using default value for parameter ' prefix params{ii} ': ' str]);
            end
        end
  end
end

% Check no options were provided that were not present in defaults
if errorwithunknownopts
    optfieldnames = fieldnames(opts);

    Idx = [];
    for ii = 1:length(optfieldnames)
        if ~isfield(defaults,optfieldnames{ii})
            Idx = [Idx ii];
        end
    end
    
    if ~isempty(Idx)
        disp('Some options provided were not present in the defaults:')
        for jj = 1:length(Idx)
            disp([prefix optfieldnames{jj}])
        end
        error('Cannot assign options - exiting!')
    end
end