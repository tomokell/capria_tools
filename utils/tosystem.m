% Calls a system command after first echoing it and printing the result to
% the screen and also printing any resulting messages
%
% Tom Okell, June 2022
%
% tosystem(cmd)
function tosystem(cmd)

% Echo the command (so wildcards etc. are replaced) then display the result
[~,cmd_msg] = system(['echo ' cmd]);
disp(cmd_msg);

% Run using builtin which forces output to be displayed
builtin('system',cmd);