%% TASK B - Ionospheric TEC: observation vs global model
% Files: ionospheredata1.dat, ionospheredata2.dat  (Time, Obs, Model) [TECU]
clear; clc; close all;

files = {'ionospheredata1.dat','ionospheredata2.dat'};
res   = struct();

for k = 1:2
    % ---------- read ----------
    fid = fopen(files{k},'r');
    fgetl(fid);                                   % skip header line
    C = textscan(fid,'%s %s %f %f');
    fclose(fid);
    t   = datetime(strcat(C{1},{' '},C{2}),'InputFormat','M/d/yyyy HH:mm:ss');
    obs = C{3};  mod = C{4};
    n   = numel(obs);

    % ---------- B1. time series of the whole span ----------
    figure('Name',sprintf('Time series data %d',k));
    plot(t,obs,'b-','LineWidth',1); hold on
    plot(t,mod,'r-','LineWidth',1); grid on
    xlabel('Time [UTC]'); ylabel('TEC [TECU]');
    legend('Observed TEC','Global model TEC','Location','best');
    title(sprintf('TEC time series - %s',files{k}));

    % ---------- B2. scatter plot + correlation ----------
    R = corrcoef(mod,obs); R = R(1,2);

    % ---------- B3. linear model  obs = a*model + b  (least squares) ----------
    x = mod; y = obs;
    A = [x, ones(n,1)];
    p = A\y;   a = p(1); b = p(2);
    yfit = A*p;
    e    = y - yfit;
    SSE  = sum(e.^2);                              % sum of squared errors
    dof  = n-2;
    SEE  = sqrt(SSE/dof);                          % standard error of estimate
    Sxx  = sum((x-mean(x)).^2);
    se_a = SEE/sqrt(Sxx);                          % std error of slope
    se_b = SEE*sqrt(1/n + mean(x)^2/Sxx);          % std error of intercept
    R2   = 1 - SSE/sum((y-mean(y)).^2);

    % differences (observation - model, no regression)
    d      = obs - mod;
    d_mean = mean(d);  d_std = std(d);
    % differences after regression = residuals e
    RMSE   = sqrt(mean(d.^2));

    figure('Name',sprintf('Scatter data %d',k));
    scatter(x,y,6,'filled','MarkerFaceAlpha',0.4); hold on
    xl = [min([x;y]) max([x;y])];
    plot(xl,xl,'k--','LineWidth',1.2);                         % 1:1 line
    plot(xl,a*xl+b,'r-','LineWidth',1.5); grid on; axis equal
    xlabel('Global model TEC [TECU]'); ylabel('Observed TEC [TECU]');
    legend('Data','1:1 line',sprintf('y = %.3fx + %.3f',a,b),'Location','best');
    title(sprintf('%s   R = %.3f,  R^2 = %.3f',files{k},R,R2));

    figure('Name',sprintf('Residuals data %d',k));
    subplot(1,2,1); plot(t,e,'k.'); grid on
    xlabel('Time [UTC]'); ylabel('Residual [TECU]'); title('Regression residuals');
    subplot(1,2,2); histogram(d,40); grid on
    xlabel('Obs - Model [TECU]'); ylabel('Count'); title('Differences (Obs - Model)');

    % ---------- print ----------
    fprintf('\n===== %s (n = %d) =====\n',files{k},n);
    fprintf('Correlation R              = %.4f   (R^2 = %.4f)\n',R,R2);
    fprintf('Model: obs = a*model + b\n');
    fprintf('  a = %.4f  +/- %.4f (std err)\n',a,se_a);
    fprintf('  b = %.4f  +/- %.4f (std err)\n',b,se_b);
    fprintf('SSE                        = %.2f TECU^2\n',SSE);
    fprintf('SEE = sqrt(SSE/(n-2))      = %.4f TECU\n',SEE);
    fprintf('Mean of differences (O-M)  = %.4f TECU\n',d_mean);
    fprintf('Std of differences (O-M)   = %.4f TECU\n',d_std);
    fprintf('RMSE (O-M)                 = %.4f TECU\n',RMSE);

    res(k).file=files{k}; res(k).R=R; res(k).R2=R2; res(k).a=a; res(k).b=b;
    res(k).se_a=se_a; res(k).se_b=se_b; res(k).SSE=SSE; res(k).SEE=SEE;
    res(k).dmean=d_mean; res(k).dstd=d_std; res(k).RMSE=RMSE;
end

%% Summary table
Tsum = struct2table(res);
disp(Tsum);

%% B4. Discussion hints (interpret with the figures/table above)
fprintf(['\n=== DISCUSSION GUIDE ===\n' ...
 '- Time series: look for diurnal cycle (TEC peaks around local afternoon, minimum near\n' ...
 '  pre-dawn), day-to-day variability, and the model being a smooth curve while observations\n' ...
 '  show sudden jumps/dips (storm / irregularities / data discontinuities).\n' ...
 '- R close to 1 : observation and model vary together (same pattern); slope a != 1 or b != 0\n' ...
 '  means scale / offset bias; mean of differences quantifies systematic bias.\n' ...
 '- Small SEE and std of differences = good agreement; large values = model cannot\n' ...
 '  reproduce the observed variability (e.g. phase shift in time or local effects).\n' ...
 '- Compare datasets: the one with higher R, a closer to 1, b closer to 0 and smaller SEE is\n' ...
 '  better represented by the global model. Standard errors of a,b show coefficient reliability.\n' ...
 '- Note: samples are autocorrelated (30 s sampling), so the standard errors are optimistic.\n']);
