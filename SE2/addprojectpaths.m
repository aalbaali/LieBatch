function addprojectpaths()
    % Add project paths.
    
    % Add paths from the root directory
    run ../addprojectpaths.m;
    
    % Add paths to the Initialize repo
    addpath( 'Initialization\');
end