function [ParmFit, X_sample, Y_sample] = mixedwblplot(x)
%WBLPLOT Weibull probability plot.
%   H = WBLPLOT(X) displays a Weibull probability plot of the data in X. For
%   matrix, X, WBLPLOT displays a plot for each column. H is a handle to the
%   plotted lines.
%   
%   The purpose of a Weibull probability plot is to graphically assess whether
%   the data in X could come from a Weibull distribution. If the data are
%   Weibull the plot will be linear. Other distribution types will introduce
%   curvature in the plot. WBLPLOT uses midpoint probability plotting
%   positions. Use PROBPLOT when the data included censored observations.
%
%   Use the data cursor to read precise values, observation numbers, and  
%   the value of the observations projected on to the reference line. 
%
%   See also PROBPLOT, NORMPLOT.

%   Copyright 1993-2010 The MathWorks, Inc. 


if size(x,1)==1
    x = x';
end
[n, m] = size(x);

[sx, originds]= sort(x);
minx  = min(sx(:));
maxx  = max(sx(:));
if isnan(minx) % Data all NaNs, setting arbitrary limits.
    minx = 1;
    maxx = 10;
end
range = maxx-minx;

if range > 0
  minxaxis  = 0;
  maxxaxis  = maxx+0.025*range;
else
  minxaxis  = minx - 1;
  maxxaxis  = maxx + 1;
end

% Use the same Y vector if all columns have the same count
if (~any(isnan(x(:))))
   eprob = ((1:n)' - 0.5)./n;
else
   nvec = sum(~isnan(x));
   eprob = repmat((1:n)', 1, m);
   eprob = (eprob-.5) ./ repmat(nvec, n, 1);
   eprob(isnan(sx)) = NaN;
   n = max(nvec);  % sample size for setting axis limits
end
y  = log(log(1./(1-eprob)));
if (size(y,2) < m)
   y = y(:, ones(1,m));
end
if n>0
    minyaxis  = log(log(1./(1 - 0.25 ./n)));
    maxyaxis  = log(log(1./(0.25 ./n)));
else
    minyaxis = -4; % Data all NaNs, setting arbitrary limits.
    maxyaxis = 1;
end


p     = [0.001 0.003 0.01 0.02 0.05 0.10 0.25 0.5...
         0.75 0.90 0.96 0.99 0.999];

label = {'0.001','0.003', '0.01','0.02','0.05','0.10','0.25','0.50', ...
         '0.75','0.90','0.96','0.99', '0.999'};

tick  = log(log(1./(1-p)));

q1x = prctile(x,25);
q3x = prctile(x,75);
q1y = prctile(y,25);
q3y = prctile(y,75);

qx = [q1x; q3x];
qy = [q1y; q3y];

% Fit the mixed Weibull distribution: base, transition, and tail parts
ReLUfcn = @(x) max(real(x),0);
model = @(P,x) P(1) + P(2)*ReLUfcn(log(x)) + (P(3)-P(2))*ReLUfcn(log(x/P(4)))+(P(5)-P(3))*ReLUfcn(log(x/P(6)));
modelseg1 = @(P,x) P(1) + P(2)*real(log(x));
modelseg2 = @(P,x) P(1) + P(2)*real(log(P(4))) + P(3)*real(log(x/P(4)));
modelseg3 = @(P,x) P(1) + P(2)*real(log(P(4))) + P(3)*real(log(P(6)/P(4))) + P(5)*real(log(x/P(6)));
P0 = [-85, 50, -20, prctile(x,82), -30, prctile(x,95)];
ParmFit = lsqcurvefit(model,P0,sx,y);
K = [ParmFit(2), ParmFit(3), ParmFit(5)]
Lambda = exp(-[modelseg1(ParmFit,1) modelseg2(ParmFit,1) modelseg3(ParmFit,1)]./K)
Mean = Lambda.*gamma(1+1./K)
% modelpred = model(ParmFit,sx);

b = zeros(m,2);
mx = [max(minx*0.95, minx-10) min(ParmFit(4)*1.1, ParmFit(4)+50) ...
        max(ParmFit(4)*0.93,ParmFit(4)-30) min(ParmFit(6)*1.1,ParmFit(6)+50) ...
        max(ParmFit(6)*0.9,ParmFit(6)-50) min(maxx*1.1, maxx+50)];
my = zeros(m,2);

my1 = modelseg1(ParmFit, mx(1:2));
my2 = modelseg2(ParmFit, mx(3:4));
my3 = modelseg3(ParmFit, mx(5:6));

% Calculate reliability estimation table
EPS=10.^(-[1:10])
Estimation=ParmFit(6)*exp((log(log(1./EPS))-model(ParmFit,ParmFit(6)))./ParmFit(5))

X_sample = 0:0.01:(max(x)+5*var(x));
Y_sample = -1./(exp(exp(model(ParmFit,X_sample))))+1./(exp(exp(model(ParmFit,X_sample-0.01))));
Y_sample(isnan(Y_sample))=0;

% Plot data and corresponding reference lines in the same color,
% following the default color order.  Plot reference line first, 
% followed by the data, so that data will be on top of reference line.
%hrefends = line(mx,my,'LineStyle','-.','Marker','none');
%hrefmid = line(qx,qy,'LineStyle','-','Marker','none');
% href1 = line(mx(1:2), my1,'LineStyle','-','Marker','none','Color','r');
% href2 = line(mx(3:4), my2,'LineStyle','-','Marker','none','Color','r');
 href3 = line(mx(5:6), my3,'LineStyle','-','Marker','none','Color','r');
%hdat = line(sx(IsIdle(originds)),y(IsIdle(originds)),'LineStyle','none','Marker','+','Color','b');
hdat = line(sx,y,'LineStyle','none','Marker','+','Color','b');
%hdatbusy = line(sx(~IsIdle(originds)),y(~IsIdle(originds)),'LineStyle','none','Marker','+','Color','m');

if m==1
    set(hdat,'MarkerEdgeColor','b');
%    set([hrefmid, hrefends],'Color','r');
end
if nargout>0
%    h = [hdat;hrefmid;hrefends];
end

if m>1
for i=1:m
    % Set custom data cursor on data
    hB = hggetbehavior(hdat(i),'datacursor');
    set(hB,'UpdateFcn',{@wblplot_datatip_callback,b(i,:)});
    % Disable datacursor on reference lines
    %hB = hggetbehavior(hrefends(i),'datacursor');
    %set(hB,'Enable',false);
    %hB = hggetbehavior(hrefmid(i),'datacursor');
    %set(hB,'Enable',false);
    if m>1
        setappdata(hdat(i),'group',i);
    end
    if ~isempty(originds)
        setappdata(hdat(i),'originds',originds(:,i));
    end
end
end

set(gca,'YTick',tick,'YTickLabel',label,'XScale','log');
set(gca,'YLim',[minyaxis maxyaxis],'XLim',[minxaxis maxxaxis]);
%xlabel(getString(message('stats:probplot:Data')));
xlabel('ET / us');
ylabel(getString(message('stats:probplot:Probability')));
title(getString(message('stats:probplot:WeibullProbPlot')));

grid on; hold on;
%plot(sx,modelpred,'r-','LineWidth',1.5); hold off;



% ----------------------------
function datatip_txt = wblplot_datatip_callback(obj,evt,b)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

x = pos(1);
y = pos(2);

% Compute position to display.
yper = 1-1./exp(exp(y));
yperexp = 1-1./exp(exp(polyval(b,log(x))));
xexp = exp(polyval([1/b(1),-b(2)/b(1)],y));

% Get the original row number of the selected point.
originds = getappdata(target,'originds');
origind = originds(ind);

% Get the group number, which is set if more than one.
group = getappdata(target,'group');
 
% Generate text to display.
datatip_txt = {
    sprintf('%s: %s',getString(message('stats:probplot:Data')),num2str(x)),...
    sprintf('%s: %s',getString(message('stats:probplot:Probability')),num2str(yper)),...
    ''
    };
datatip_txt{end+1} = sprintf('%s: %s',...
    getString(message('stats:probplot:Observation')),num2str(origind));
if ~isempty(group)
    datatip_txt{end+1} = sprintf('%s: %s',...
        getString(message('stats:probplot:Group')),num2str(group));
end

datatip_txt{end+1} = '';
datatip_txt{end+1} = sprintf('%s %s: %s',...
    getString(message('stats:probplot:ReferenceLine')),...
    getString(message('stats:probplot:Data')),...
    num2str(xexp));
datatip_txt{end+1} = sprintf('%s %s: %s',...
    getString(message('stats:probplot:ReferenceLine')),...
    getString(message('stats:probplot:Probability')),...
    num2str(yperexp));
