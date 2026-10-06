%% TASK 1 - Time series plots and outlier detection
% Datasets: sealevelobs.dat (hourly sea level, Oslo, cm)
%           updeformation.dat (daily GNSS CORS CTGL Up component, m)
% Put the .dat files in the same folder as this script.
clear; clc; close all;

%% ===================== 1. SEA LEVEL =====================
fid = fopen('sealevelobs.dat','r');
C   = textscan(fid,'%s %f','CommentStyle','#');
fclose(fid);
tstr = C{1};  sl = C{2};
t_sl = datetime(cellfun(@(s) s(1:19),tstr,'UniformOutput',false), ...
                'InputFormat','yyyy-MM-dd''T''HH:mm:ss');

% --- (1) plot raw series
figure('Name','Sea level raw');
plot(t_sl, sl, 'b.-'); grid on
xlabel('Time (local, UTC+2)'); ylabel('Sea level variation [cm]');
title('Oslo water level observations (raw)');

% --- (2) remove deterministic part: mean + tidal harmonics (least squares)
th = hours(t_sl - t_sl(1));                       % time in hours
per = [12.4206 12.0000 12.6583 11.9672 ...        % M2 S2 N2 K2
       23.9345 25.8193 24.0659 26.8684 ...        % K1 O1 P1 Q1
       327.86 661.31];                            % Mf (13.66 d), Mm (27.55 d)
A = ones(numel(th),1);
for k = 1:numel(per)
    A = [A, cos(2*pi*th/per(k)), sin(2*pi*th/per(k))]; %#ok<AGROW>
end
beta   = A \ sl;
sl_det = A*beta;
res_sl = sl - sl_det;                             % residual (stochastic + outliers)

figure('Name','Sea level deterministic');
subplot(2,1,1)
plot(t_sl, sl,'b.-', t_sl, sl_det,'r-'); grid on
xlabel('Time'); ylabel('Sea level [cm]'); legend('Observed','Harmonic tide model');
subplot(2,1,2)
plot(t_sl, res_sl,'k.-'); grid on
xlabel('Time'); ylabel('Residual [cm]'); title('Residual after removing tide');

% --- outlier detection on residuals
[out2s_sl, lim2s_sl] = twosigma(res_sl);
[outIQR_sl, limIQR_sl] = iqrrule(res_sl);
fprintf('\n=== SEA LEVEL ===\n');
fprintf('N = %d | residual mean = %.2f cm, std = %.2f cm\n', numel(res_sl), mean(res_sl), std(res_sl));
fprintf('Two-sigma rule : %d outliers (limits %.2f .. %.2f cm)\n', sum(out2s_sl), lim2s_sl);
fprintf('IQR rule       : %d outliers (limits %.2f .. %.2f cm)\n', sum(outIQR_sl), limIQR_sl);
fprintf('Outlier times (IQR):\n'); disp(t_sl(outIQR_sl)'); disp(sl(outIQR_sl)');

% skewness/kurtosis help decide the method
fprintf('Skewness = %.2f, kurtosis = %.2f (normal: 0 and 3)\n', skew(res_sl), kurt(res_sl));

figure('Name','Sea level outliers');
subplot(2,1,1)
plot(t_sl,res_sl,'k.'); hold on
plot(t_sl(outIQR_sl),res_sl(outIQR_sl),'ro','MarkerSize',8);
plot(t_sl(out2s_sl),res_sl(out2s_sl),'gx','MarkerSize',8);
yline(limIQR_sl(1),'r--'); yline(limIQR_sl(2),'r--');
yline(lim2s_sl(1),'g--');  yline(lim2s_sl(2),'g--'); grid on
xlabel('Time'); ylabel('Residual [cm]');
legend('Residual','IQR outliers','2\sigma outliers','IQR limits','','2\sigma limits');
title('Outlier detection on sea-level residuals');
subplot(2,1,2)
histogram(res_sl,40); hold on; xline(limIQR_sl(1),'r--'); xline(limIQR_sl(2),'r--'); xline(lim2s_sl(1),'g--'); xline(lim2s_sl(2),'g--'); grid on; xlabel('Residual [cm]'); ylabel('Count'); title('Histogram of residuals (red = IQR, green = 2\sigma limits)');

% --- remove outliers (IQR rule chosen, see discussion) and replot
sl_clean = sl(~outIQR_sl);  t_sl_clean = t_sl(~outIQR_sl);
figure('Name','Sea level cleaned');
plot(t_sl_clean, sl_clean,'b.-'); grid on
xlabel('Time (local, UTC+2)'); ylabel('Sea level variation [cm]');
title('Oslo water level - outliers removed');

%% ===================== 2. GNSS UP COMPONENT =====================
T = readtable('updeformation.dat','Delimiter',',','ReadVariableNames',true);
t_up = datetime(T{:,1},'InputFormat','yyyy-MM-dd HH:mm:ss');
up   = T{:,2};

% --- (1) plot raw
figure('Name','Up raw');
plot(t_up, up,'b.'); grid on
xlabel('Time [UTC]'); ylabel('Relative Up displacement [m]');
title('CTGL GNSS CORS Tegal - Up component (raw)');

% --- (2) deterministic part: cubic polynomial + annual + semi-annual
td = days(t_up - t_up(1));  tn = td/max(td);       % normalised for conditioning
B  = [ones(size(tn)) tn tn.^2 tn.^3 ...
      cos(2*pi*td/365.25) sin(2*pi*td/365.25) ...
      cos(4*pi*td/365.25) sin(4*pi*td/365.25)];
% robust fit (iteratively reweighted) so that outliers don't bias the model
w = ones(size(up));
for it = 1:10
    W  = diag(sparse(w));
    bb = (B'*W*B)\(B'*W*up);
    r  = up - B*bb;
    s  = 1.4826*median(abs(r-median(r)));          % MAD scale
    w  = 1./max(1,abs(r)/(3*s));                   % Huber-like weights
end
up_det = B*bb;
res_up = up - up_det;

figure('Name','Up deterministic');
subplot(2,1,1)
plot(t_up,up,'b.',t_up,up_det,'r-','LineWidth',1.2); grid on
xlabel('Time [UTC]'); ylabel('Up [m]'); legend('Observed','Trend + seasonal model');
subplot(2,1,2)
plot(t_up,res_up,'k.'); grid on
xlabel('Time [UTC]'); ylabel('Residual [m]'); title('Residual after removing trend & seasonal');

[out2s_up, lim2s_up]   = twosigma(res_up);
[outIQR_up, limIQR_up] = iqrrule(res_up);
fprintf('\n=== GNSS UP ===\n');
fprintf('N = %d | residual mean = %.4f m, std = %.4f m\n', numel(res_up), mean(res_up), std(res_up));
fprintf('Two-sigma rule : %d outliers (limits %.4f .. %.4f m)\n', sum(out2s_up), lim2s_up);
fprintf('IQR rule       : %d outliers (limits %.4f .. %.4f m)\n', sum(outIQR_up), limIQR_up);
fprintf('Skewness = %.2f, kurtosis = %.2f\n', skew(res_up), kurt(res_up));
fprintf('Largest IQR outliers:\n');
[~,ix] = sort(abs(res_up),'descend'); ix = ix(1:8);
disp(table(t_up(ix), up(ix), res_up(ix), outIQR_up(ix), 'VariableNames',{'Time','Up_m','Res_m','IsOutlier_IQR'}));

figure('Name','Up outliers');
subplot(2,1,1)
plot(t_up,res_up,'k.'); hold on
plot(t_up(outIQR_up),res_up(outIQR_up),'ro','MarkerSize',8);
plot(t_up(out2s_up),res_up(out2s_up),'gx','MarkerSize',8);
yline(limIQR_up(1),'r--'); yline(limIQR_up(2),'r--');
yline(lim2s_up(1),'g--');  yline(lim2s_up(2),'g--'); grid on
xlabel('Time [UTC]'); ylabel('Residual [m]');
legend('Residual','IQR outliers','2\sigma outliers');
title('Outlier detection on GNSS Up residuals');
subplot(2,1,2)
histogram(res_up,40); hold on; xline(limIQR_up(1),'r--'); xline(limIQR_up(2),'r--'); xline(lim2s_up(1),'g--'); xline(lim2s_up(2),'g--'); grid on; xlabel('Residual [m]'); ylabel('Count'); title('Histogram of residuals (red = IQR, green = 2\sigma limits)');

% --- remove outliers & replot
up_clean = up(~outIQR_up); t_up_clean = t_up(~outIQR_up);
figure('Name','Up cleaned');
plot(t_up_clean, up_clean,'b.'); grid on
xlabel('Time [UTC]'); ylabel('Relative Up displacement [m]');
title('CTGL GNSS Up - outliers removed');

%% ===================== 3. DISCUSSION (printed) =====================
fprintf('\n=== DISCUSSION ===\n');
fprintf(['Two-sigma rule assumes roughly Gaussian residuals and flags ~5%% of normal data by\n' ...
         'construction; mean/std are themselves inflated by outliers (masking).\n' ...
         'IQR rule uses quartiles (robust) and makes no normality assumption, so it is better\n' ...
         'for data with heavy tails / isolated spikes (GNSS Up: spikes of 0.1 m vs noise of ~0.01 m).\n' ...
         'For sea level, the residual after tide removal is dominated by meteorological (surge)\n' ...
         'signal and is autocorrelated; compare the counts above and check skewness/kurtosis.\n']);

%% ===================== local functions =====================
function [isOut, lims] = twosigma(x)
    m = mean(x); s = std(x);
    lims = [m-2*s, m+2*s];
    isOut = x < lims(1) | x > lims(2);
end
function [isOut, lims] = iqrrule(x)
    q = prctile(x,[25 75]); i = q(2)-q(1);
    lims = [q(1)-1.5*i, q(2)+1.5*i];
    isOut = x < lims(1) | x > lims(2);
end
function s = skew(x), x = x(:); s = mean((x-mean(x)).^3)/std(x,1)^3; end
function k = kurt(x), x = x(:); k = mean((x-mean(x)).^4)/std(x,1)^4; end
