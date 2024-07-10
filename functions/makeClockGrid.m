% Modified from: Perfect Polar Plots 'makePolarGrid'
% https://www.mathworks.com/matlabcentral/fileexchange/73967-perfect-polar-plots
%--------------------------------------------------------------------------


function ph = makeClockGrid(ATicks, RTicks, time_type, varargin)

% default clock params
t_min= 0;  t_max= 24; 
ADir=     'clockwise';  % 'clockwise'/'ccw' or 'counterclockwise'/'ccw'
ALim=     [0, 360];     % span of circle grid
AZeroPos= 'top';        % initial 0deg position 'top',left','bottom','right'
AUnits=   '';           % (optional) ignored for some, used only for spokeLabelSchemes
AFontSize= 10;          % (optional) change font size


RLablScheme=  'normal'; % various schemes available
RLim=         []; 
RMinorTicks=  [];       % (optional) extra circles, but won't have labels
RLablAng=     15;       % Angle where the radial tick labels will go
RLabelFormat= '%.1f';   % Format for radial labels
RFontSize=    10;       % (optional) change font size
RUnits=       '';       % (optional) will be applied to last radial label


% clock plot properties
rotateGrid = 0;  % Rotate grid arbitrary amount (degrees, positive in direction of ADir)

majorGridLineStyle  = '-';             % major grid line style
majorGridLineWidth  = 0.25;            % major grid LineWidth
majorGridColor      = [0.5, 0.5, 0.5]; % major grid color

minorGridLineStyle  = ':';             % minor grid line style
minorGridLineWdith  = 0.1;             % minor grid LineWidth
minorGridColor      = [0.5, 0.5, 0.5]; % minor grid color

outerBoxLineStyle   = '-';             % outer bounding box line style
outerBoxLineWidth   = 1;               % outer bounding box LineWidth
outerBoxColor       = 'k';             % outer bounding box color


% parse variable inputs
Aminor= strcmpi(varargin, 'AMinorTicks'); 
if any(Aminor)
    AMinorTicks= varargin{find(Aminor) + 1};
else
    AMinorTicks= 5;     % (optional) extra spokes, but won't have labels
end

minorGrid= strcmpi(varargin, 'minorGrid');
if any(minorGrid)
    showMinorGrid= varargin{find(Aminor) + 1};
    switch showMinorGrid
        case 'off'
            showMinorGrid= 0;
        otherwise
            showMinorGrid= 1;   % minor grid on by default
    end
else
    showMinorGrid= 1;   % minor grid on by default
end

nSubplot= strcmpi(varargin, 'nSubplots');
if any(nSubplot)
    txt_scale= varargin{find(nSubplot) + 1};
else
    txt_scale= 1;     % (optional) extra spokes, but won't have labels
end


% If a single value is given for ATicks, then its the increment.
if length(ATicks) == 1
  ATicks = ALim(1) : ATicks : ALim(2);
end
ATicks = ATicks(:)';  % Make row vector

n_hrs= length(ATicks)-1; 

% make clock labels
labels= makeClockLabels(t_min, t_max, n_hrs, time_type);

% make row vector and reduce overlap with Major ATicks
AMinorTicks= ALim(1) : AMinorTicks : ALim(2);
AMinorTicks= AMinorTicks(~ismember(AMinorTicks, ATicks));


if isempty(RLim)
  RLim = [RTicks(1), RTicks(end)];
end

axATicks      = convStdPolarAngle(ATicks, ADir, AZeroPos, rotateGrid);
axAMinorTicks = convStdPolarAngle(AMinorTicks, ADir, AZeroPos, rotateGrid);
axRLabelAngle = convStdPolarAngle(RLablAng, ADir, AZeroPos, rotateGrid);

ph.ALim          = ALim;
ph.ATicks        = ATicks;
ph.AFontSize     = AFontSize; 
ph.axATicks      = axATicks;
ph.AMinorTicks   = AMinorTicks;
ph.axAMinorTicks = axAMinorTicks;
ph.axRLabelAngle = axRLabelAngle;
ph.RTicks        = RTicks;
ph.RLim          = RLim;
ph.RMinorTicks   = RMinorTicks;
ph.RLabelFormat  = RLabelFormat;
ph.RFontSize     = RFontSize; 
ph.RUnits        = RUnits; 



%% Set Clock Figure Properties --------------------------------------------

hold on
cax = gca;
set(gcf,'color','w');
%set(gca,'Color','none');
% set(gcf,'Position',[200 50 900 700]); % this size is abitrary
view(cax, 2);  % set view to 2-D

set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
set(get(cax, 'XLabel'), 'Visible', 'on');
set(get(cax, 'YLabel'), 'Visible', 'on');

% Set size of plot area 
% This may need to be adjusted if doesnt allow subplots, but author finds 
% it works fine in most cases.
axis(cax, RLim(2) * [-1, 1, -1.15, 1.15]);

ph.figHnd = gcf;
ph.axHnd  = cax;

% Define a circle. These are used to help draw
% the various grid lines etc.
if ALim(1) > ALim(2)
  th = ALim(1) : 1 : (ALim(2)+360);
else
  th = ALim(1) : 1 : ALim(2);
end
th = convStdPolarAngle(th, ADir, AZeroPos, rotateGrid) * pi/180;  
xunit = cos(th);
yunit = sin(th);


%% Draw Minor Spoke Grid Lines --------------------------------------------

if showMinorGrid && ~isempty(AMinorTicks)
  px = [RLim(1)*cosd(axAMinorTicks)', RLim(2)*cosd(axAMinorTicks)']';
  py = [RLim(1)*sind(axAMinorTicks)', RLim(2)*sind(axAMinorTicks)']';
  
  hh = line(px,py,...
    'LineStyle', minorGridLineStyle,...
    'Color',     minorGridColor,...
    'LineWidth', minorGridLineWdith, ...
    'HandleVisibility', 'off',...
    'Parent', cax);
  
  ph.AMinorGrid = hh;   clear px py hh
else
  ph.AMinorGrid = [];
end


%% Draw Major Spoke Grid Lines ---------------------------------------------

if ~isempty(ATicks)
  dax= axATicks;
  
  px= [RLim(1)*cosd(dax)', RLim(2)*cosd(dax)']';
  py= [RLim(1)*sind(dax)', RLim(2)*sind(dax)']';
  
  hh = line(px, py, ...
    'LineStyle', majorGridLineStyle,...
    'Color',     majorGridColor,...
    'LineWidth', majorGridLineWidth, ...
    'HandleVisibility', 'off',...
    'Parent', cax);
  
  ph.AMajorGrid = hh;   clear px py dax
else
  ph.AMajorGrid = [];
end


%% Draw Minor Circle Grid Lines -------------------------------------------

if showMinorGrid && ~isempty(RMinorTicks)
  [xm, ym] = meshgrid(xunit, RMinorTicks);
  px= xm .* ym;
  clear xm ym
  
  [xm, ym] = meshgrid(yunit, RMinorTicks);
  py= xm .* ym;
  clear xm ym

  hh = line(px', py', ...
      'LineStyle', minorGridLineStyle, ...
      'Color',     minorGridColor, ...
      'LineWidth', minorGridLineWdith, ...
      'HandleVisibility', 'off', ...
      'Parent', cax);
    
  ph.RMinorGrid = hh;   clear px py hh
else
  ph.RMinorGrid = [];
end


%% Draw Major Circle Grid Lines -------------------------------------------

crad = RTicks;

% Don't draw grid circle at center or outer boundaries
id = (crad == 0) | ismember(crad, RLim);
crad(id) = [];
clear id

if ~isempty(crad)
  [xm,ym] = meshgrid(xunit,crad);
  px = xm .* ym;
  clear xm ym
  
  [xm,ym] = meshgrid(yunit,crad);
  py = xm .* ym;
  clear xm ym
  
  hh = line(px', py', ...
    'LineStyle', majorGridLineStyle, ...
    'Color',     majorGridColor, ...
    'LineWidth', majorGridLineWidth, ...
    'HandleVisibility', 'off', ...
    'Parent', cax);
  
  ph.RMajorGrid = hh;   clear px py hh
else
  ph.RMajorGrid = [];
end


%% Draw Outer Boundary-----------------------------------------------------

% Draw outer circle
ph.OuterBox = line(xunit * RLim(2), yunit * RLim(2), ...
  'LineStyle', outerBoxLineStyle, ...
  'Color',     outerBoxColor, ...
  'LineWidth', outerBoxLineWidth, ...
  'HandleVisibility', 'off', ...
  'Parent', cax);

% Draw inner circle
if RLim(1) > 0
  
  ph.OuterBox(end+1,1) = line(xunit * RLim(1), yunit * RLim(1), ...
    'LineStyle', outerBoxLineStyle, ...
    'Color',     outerBoxColor, ...
    'LineWidth', outerBoxLineWidth, ...
    'HandleVisibility', 'off', ...
    'Parent', cax);
  
end


%% Label Spoke Ticks ------------------------------------------------------

for k= 1:length(ATicks)-1
  tickAngle = ATicks(k);
  axAngle   = axATicks(k);
  tmplabel  = labels{k}; 
  
  if axAngle == axRLabelAngle
    continue
  end

  ha = annotateNormal(tickAngle, axAngle, RLim(2), tmplabel, txt_scale, ph);
  
  ph.ALabels(k) = ha;      clear tickAngle axAngle ha
end


%% Label Radial Ticks ------------------------------------------------------

for r = 1:length(RTicks)
  thisRad = RTicks(r);
  
  switch RLablScheme
    case 'none'
      hr = [];    
    case 'normal'
      hr = getRLabel_Normal(thisRad, txt_scale, ph);
    case 'outer'  
      hr = getRLabel_Outer(thisRad, txt_scale, ph);
    otherwise
      error(['unrecognized RLabelScheme: ',RLablScheme]);
  end
      
  if ~isempty(hr)
    ph.RLabels(r) = hr;
  end
  clear radStr hr
  
end

fn = fieldnames(ph);
[~,id] = sort(lower(fn));
ph = orderfields(ph,fn(id));
clear fn id


end  % END MAIN FUNCTION



%% LOCAL FUNCTIONS


%% function shift2radius --------------------------------------------------
function labHnd = shift2radius(labHnd, rclose, axAngle, angShiftFlag)
% shifts a text box so it the closest its extent box is rt from a circle of
% radius rt from the origin

if nargin<4
  angShiftFlag = false;
end

rt = rclose;

for k = 1:5

b = labHnd.Extent;

p = [b(1), b(2);...
  b(1)+b(3)/2, b(2);...
  b(1)+b(3), b(2);...
  b(1)+b(3), b(2)+b(4)/2;...
  b(1)+b(3), b(2)+b(4);...
  b(1)+b(3)/2, b(2)+b(4);...
  b(1), b(2)+b(4);...
  b(1), b(2)+b(4)/2];

% Just checking these match
%pang = wrap360(atan2d(mean(p(:,2)),mean(p(:,1))));

[d,id] = min(sqrt(p(:,1).^2 + p(:,2).^2));

if k==1 && angShiftFlag
dtheta = wrap360(atan2d(p(id,2),p(id,1))) - wrap360(axAngle);
if abs(dtheta)>3 && abs(dtheta)<10
  axAngle = axAngle - dtheta/2;
end
end

rt = rclose*(rt/d);
labHnd.Position = [rt*cosd(axAngle), rt*sind(axAngle)];
  
%chk(k) = axTheta;

end

end % shift2radius


%% function getRLabel_Normal ----------------------------------------------
function hr = getRLabel_Normal(thisRad, sc, ph)

% dr is tiny shift to put the labels just outside of each grid circle
hr = [];

if thisRad==0
  return;
end

RLim = ph.RLim;

dr = (1.005 - (.005 * (1 - 1/sc))) * RLim(2) - RLim(2);

if (thisRad == RLim(2)) && ~isempty(ph.RUnits)
  radStr = [sprintf(ph.RLabelFormat,thisRad),' ',ph.RUnits];
else
  radStr = sprintf(ph.RLabelFormat,thisRad);
end

% Positions for labels, not a slight shift
xr = (thisRad) * cosd(ph.axRLabelAngle);
yr = (thisRad) * sind(ph.axRLabelAngle);

hr = text(xr, yr, radStr, ...
  'VerticalAlignment',   'middle', ...
  'HorizontalAlignment', 'center', ...
  'HandleVisibility',    'off', ...
  'FontSize', ph.RFontSize, ...
  'Parent', ph.axHnd);

if thisRad == RLim(2)
  hr = shift2radius(hr, 1.01*RLim(2),  ph.axRLabelAngle, true);
else
  hr = shift2radius(hr, thisRad+dr, ph.axRLabelAngle, false);
end

end   % getRLabel_Normal


%% function annotateNormal -------------------------------------------------
function spkHnd = annotateNormal(~, d, rmax, clabl, sc, ph)

% Default distance from center for label
rt = (1.005 - (.005 * (1 - 1/sc))) * rmax;
labStr = clabl;
% labStr = sprintf('%.0f%s', clabl, char(176));
  
spkHnd = text(rt*cosd(d), rt*sind(d), labStr,...
  'HorizontalAlignment', 'center', ...
  'VerticalAlignment', 'middle', ...
  'HandleVisibility', 'off',...
  'FontSize', ph.AFontSize,...
  'Parent', ph.axHnd);

spkHnd = shift2radius(spkHnd, rt, d);


end % annotateSpokes1


%% function getRLabel_Outer -----------------------------------------------
function hr = getRLabel_Outer(thisRad, sc, ph)

% dr is tiny shift to put the labels just outside of each grid circle
hr = [];

if thisRad==0 && ph.fullCircle
  return
end

RLim = ph.RLim;
dr = (1.004 - (.005 * (1 - 1/sc))) * RLim(2) - RLim(2);
radStr = sprintf(ph.RLabelFormat,thisRad);

% Positions for labels, not a slight shift
xr = (thisRad) * cosd(ph.axRLabelAngle);
yr = (thisRad) * sind(ph.axRLabelAngle);

xs = dr * cosd(90+ph.axRLabelAngle);
ys = dr * sind(90+ph.axRLabelAngle);

hr = text(xr-xs, yr-ys, radStr, ...
  'VerticalAlignment',   'middle', ...
  'HorizontalAlignment', 'center', ...
  'HandleVisibility',    'off', ...
  'Rotation', ph.axRLabelAngle, ...
  'FontSize', ph.RFontSize, ...
  'Parent', ph.axHnd);

if (thisRad == RLim(2))  && ~isempty(ph.RUnits)
  dr = 1.11*RLim(2) - RLim(2);

  xr = mean(RLim) * cosd(ph.axRLabelAngle);
  yr = mean(RLim) * sind(ph.axRLabelAngle);

  xs = dr * cosd(90+ph.axRLabelAngle);
  ys = dr * sind(90+ph.axRLabelAngle);
  
  hr = text(xr-xs, yr-ys, ph.RUnits, ...
    'VerticalAlignment',   'middle', ...
    'HorizontalAlignment', 'center', ...
    'HandleVisibility',    'off', ...
    'Rotation', ph.axRLabelAngle, ...
    'FontSize', ph.RFontSize, ...
    'Parent', ph.axHnd);
  
end

end   % getRLabel_Outer