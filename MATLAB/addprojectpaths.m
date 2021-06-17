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
    
    % Add path to RandomVariable IO
    addpath('\\wsl$\Ubuntu-20.04\home\aa\Documents\Code_base\Local\RandomVariable\MATLAB_txt_IO');
end