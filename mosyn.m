function mosyn()
    % mosyn is the main function that adds necessary directories to the
    % MATLAB path and launches the Mosyn graphical user interface (GUI).
    %

    % Get the path of the Mosyn directory.
    dir = fileparts(which('Mosyn'));

    % Add necessary directories to the MATLAB path.
    addpath(dir);
    addpath([dir filesep 'utils']);
    addpath([dir filesep 'measures']);
    addpath([dir filesep 'classes']);
    addpath([dir filesep 'external']);
    addpath([dir filesep 'graphs']);
    addpath([dir filesep 'plugins']);
    addpath([dir filesep 'gui']);

    MosynGUI()

end