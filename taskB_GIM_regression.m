%% TASK B ONLY - CAS Global Ionosphere Map (IONEX), 1 Jan 2025, 30-min, 5 x 2.5 deg grid
% Regression / correlation of TEC at Tegal (CTGL) vs a northern point (12.5N,110E)
% Outlier screening is applied only to the regression residuals.
% Put CAS0OPSFIN_20250010000_01D_30M_GIM.INX in the same folder and press Run.
clear; clc; close all;
fn   = 'CAS0OPSFIN_20250010000_01D_30M_GIM.INX';
lats = 87.5:-2.5:-87.5;   lons = -180:5:180;

%% 1. READ IONEX (TEC maps in 0.1 TECU, 9999 = missing)
L = splitlines(fileread(fn));
TEC = nan(numel(lats),numel(lons),49);  RMS = TEC;  tt = NaT(49,1);
kind = '';  m = 0;  mr = 0;  i = 1;
while i <= numel(L)
    s = L{i};
    lab = ''; if length(s) > 60, lab = strtrim(s(61:end)); end
    switch lab
        case 'START OF TEC MAP', kind = 'TEC'; m = m+1;
        case 'START OF RMS MAP', kind = 'RMS'; mr = mr+1;
        case 'EPOCH OF CURRENT MAP'
            v = sscanf(s(1:36),'%d')';
            if strcmp(kind,'TEC'), tt(m) = datetime(v(1),v(2),v(3),v(4),v(5),v(6)); end
        case 'LAT/LON1/LON2/DLON/H'
            lat = str2double(s(3:8));  r = round((87.5-lat)/2.5)+1;
            vals = [];  i = i+1;
            while numel(vals) < numel(lons)
                vals = [vals, sscanf(L{i},'%d')']; i = i+1; %#ok<AGROW>
            end
            vals = double(vals);  vals(vals==9999) = NaN;  vals = vals*0.1;
            if strcmp(kind,'TEC'), TEC(r,:,m) = vals; else, RMS(r,:,mr) = vals; end
            continue
    end
    i = i+1;
end
fprintf('Read %d TEC maps, %d RMS maps\n', m, mr);

%% 2. EXTRACT TIME SERIES (nearest grid point)
[~,iT] = min(abs(lats-(-7.5)));  [~,jT] = min(abs(lons-110));    % Tegal (CTGL)
[~,iN] = min(abs(lats-12.5));    [~,jN] = min(abs(lons-110));    % northern point
tegal = squeeze(TEC(iT,jT,:));   north = squeeze(TEC(iN,jN,:));
rmsT  = squeeze(RMS(iT,jT,:));
th = hours(tt - tt(1));

%% 3. REGRESSION  y(Tegal) = a*x(North) + b
x = north;  y = tegal;
S = regstats_simple(x,y);
[oi,li] = iqrrule(S.e);                         % outliers in regression residuals
C = regstats_simple(x(~oi), y(~oi));            % refit without outliers
print_reg('ALL DATA', S);   print_reg('OUTLIERS REMOVED', C);
fprintf('Outliers removed: %d (limits %.2f .. %.2f TECU)\n', sum(oi), li);
% does correlation survive after removing the common diurnal cycle?
Ad = ones(size(th)); for h = 1:3, Ad = [Ad, cos(2*pi*h*th/24), sin(2*pi*h*th/24)]; end
rx = x - Ad*(Ad\x);  ry = y - Ad*(Ad\y);
c = corrcoef(rx,ry);  fprintf('R after removing diurnal cycle = %.3f\n', c(1,2));

figure('Name','Part B');
subplot(1,2,1); scatter(x,y,25,'filled'); hold on; scatter(x(oi),y(oi),40,'r','filled');
xl = [min([x;y]) max([x;y])]; plot(xl,xl,'k--',xl,S.a*xl+S.b,'r-','LineWidth',1.3); grid on
xlabel('GIM TEC 12.5N,110E [TECU]'); ylabel('GIM TEC Tegal [TECU]');
title(sprintf('R = %.3f, R^2 = %.3f',S.R,S.R2)); legend('Data','Outlier','1:1',sprintf('y=%.2fx+%.2f',S.a,S.b),'Location','best');
subplot(1,2,2); plot(tt,x,'b.-',tt,y,'r.-'); hold on
fill([tt;flipud(tt)],[y-rmsT;flipud(y+rmsT)],'r','FaceAlpha',.2,'EdgeColor','none'); grid on
xlabel('Time [UTC]'); ylabel('TEC [TECU]'); legend('12.5N','Tegal','Tegal \pm RMS');

figure('Name','Map'); pcolor(lons,lats,TEC(:,:,13)); shading flat; colorbar; hold on
plot(110,-7.5,'rp','MarkerSize',12,'MarkerFaceColor','r'); xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
title(['CAS GIM ' char(tt(13)) ' UTC (TECU)']);

%% ---------------- local functions ----------------
function S = regstats_simple(x,y)
    n = numel(x);  A = [x(:) ones(n,1)];  p = A\y(:);  e = y(:) - A*p;
    S.n = n; S.a = p(1); S.b = p(2); S.e = e; S.SSE = sum(e.^2); S.SEE = sqrt(S.SSE/(n-2));
    Sxx = sum((x-mean(x)).^2);  S.se_a = S.SEE/sqrt(Sxx);  S.se_b = S.SEE*sqrt(1/n + mean(x)^2/Sxx);
    r = corrcoef(x,y);  S.R = r(1,2);  S.R2 = S.R^2;
    d = y(:)-x(:);  S.dmean = mean(d);  S.dstd = std(d);
end
function print_reg(tag,S)
    fprintf('\n--- %s (n=%d) ---\n a = %.4f +/- %.4f\n b = %.4f +/- %.4f\n SSE = %.1f  SEE = %.3f TECU\n R = %.4f  R^2 = %.4f\n mean diff = %.3f  std diff = %.3f TECU\n', ...
        tag,S.n,S.a,S.se_a,S.b,S.se_b,S.SSE,S.SEE,S.R,S.R2,S.dmean,S.dstd);
end
function [isOut,lims] = iqrrule(x), q = prctile(x,[25 75]); d = q(2)-q(1); lims = [q(1)-1.5*d q(2)+1.5*d]; isOut = x<lims(1)|x>lims(2); end
