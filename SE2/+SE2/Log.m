%   Amro Al Baali
%   2-Mar-2021
function out = Log( X)
    % LOG( se2_pose) : SE2 -> se2alg coordinate computes the SE2 Log
    % function.
    
    out = se2alg.vee( SE2.logMap( X));
end