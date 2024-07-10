% Modified from: Perfect Polar Plots 'convToStdPolarAngle'
% https://www.mathworks.com/matlabcentral/fileexchange/73967-perfect-polar-plots
%--------------------------------------------------------------------------


function axAng = convStdPolarAngle(angleData, ADir, AZeroPos, rotateGrid)

% Get proper angle direction
switch lower(ADir)
  case {'clockwise', 'cw'}
    axAng = -angleData - rotateGrid;
  case {'counterclockwise', 'counter-clockware', 'ccw'}
    axAng = angleData + rotateGrid;
  otherwise
    error(['unrecognized angleDir: ', ADir]);
end

% Get proper position for zero angle position
switch lower(AZeroPos)
  case {'right', 'r'}
    axAng = axAng + 0;
  case {'top', 't'}
    axAng = axAng + 90;
  case {'left', 'l'}
    axAng = axAng + 180;
  case {'bottom', 'b'}
    axAng = axAng - 90;
  otherwise
    error(['unrecognized zeroPos: ', AZeroPos]);
end

axAng= mod(axAng, 360); 

end