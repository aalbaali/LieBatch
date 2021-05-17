function addprojectpaths()
    % Load project paths
    
    % Add my personal paths
    try 
        addmypaths();
    catch e
        warning(e.message);        
    end
    
    % Add external projects
    addpath( genpath( 'External'));
    
    % Add utils
    addpath( 'Utils');
    
    addpath( genpath( ...
        'G:\My Drive\Professional\Code_base\External\MATLAB\YAMLMatlab_0.4.3'));
end