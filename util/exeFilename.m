%%
function [exedir, exename, ext] = exeFilename()
    % Determine full path to calling file
    % In a deployed program, this will point to the application exe.
    % In regular Matlab, it will be the calling m-file.

    if isdeployed % Executable in stand-alone mode.
        [~, result] = system('path');
        exedir = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));

        args = System.Environment.GetCommandLineArgs;
        relpath = char(args(1));
        [~,exename,ext] = fileparts(relpath);

    else % MATLAB mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind+1).file;
        [exedir,exename,ext] = fileparts(relpath);
    end
end
