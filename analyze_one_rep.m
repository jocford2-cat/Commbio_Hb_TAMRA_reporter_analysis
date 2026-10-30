%%%analyze_one_rep.m%%%%

function metrics = analyze_one_rep(t, OD_no_sub, OD_with_sub, Flc_substrate_only, Flc_with_sub)
% Returns a scalar struct with fixed fields; missing/unavailable values are NaN.

% ---------- Initialize outputs (fixed schema) ----------
metrics = struct( ...
    'tau_ncc', NaN, ...
    'tau_deriv', NaN, ...
    'lead_half', NaN, ...
    'avg_rate_sg_mean', NaN, ...
    'avg_rate_sg', NaN, ...
    'avg_rate_spline', NaN, ...
    'avg_rate_spline_int', NaN, ...
    'peak_rate_sg', NaN, ...
    't_peak_center', NaN, ...
    't_peak_start', NaN, ...
    't_peak_end', NaN, ...
    'n_masked_ratio', NaN, ...
    'frac_masked_ratio', NaN, ...
    'first_valid_ratio_time', NaN ...
);

% ---------- Early checks & columnize ----------
if numel(t) < 5
    warning('analyze_one_rep: too few points'); 
    return
end
t = t(:);
OD_no_sub = OD_no_sub(:);
OD_with_sub = OD_with_sub(:);
Flc_substrate_only = Flc_substrate_only(:);
Flc_with_sub = Flc_with_sub(:);

% ---------- 1) Smooth-fit the basal OD curve (no substrate) ----------
smoothParam_OD = 0.2;
pp_baseline_OD = csaps(t, OD_no_sub, smoothParam_OD);
OD_no_sub_fit  = fnval(pp_baseline_OD, t);

% ---------- 2) Compute baseline-corrected biomass ----------
OD_corr  = OD_with_sub - OD_no_sub_fit;
resid_OD = OD_no_sub   - OD_no_sub_fit;
OD_thresh = max(3*std(resid_OD), 0);
OD_corr_s = movmean(OD_corr, 3);
OD_corr_masked = OD_corr_s;
%OD_corr_masked(OD_corr_s <= OD_thresh) = NaN;

% ---------- 3) Smooth-fit the baseline fluorescence drift (substrate only) ----------
smoothParam_Flc   = 0.7;
pp_baseline_Flc   = csaps(t, Flc_substrate_only, smoothParam_Flc);
Flc_substrate_fit = fnval(pp_baseline_Flc, t);

% ---------- 4) Compute baseline-corrected fluorescence ----------
Flc_corr  = Flc_with_sub - Flc_substrate_fit;
resid_Flc = Flc_substrate_only - Flc_substrate_fit;
Flc_thresh = max(3*std(resid_Flc), 0);
Flc_corr_s = movmean(Flc_corr, 3);
Flc_corr_masked = Flc_corr_s;
%Flc_corr_masked(Flc_corr_s <= Flc_thresh) = NaN;

% ---------- 5) Ratio ----------
ratio_Flc_ODcorr = Flc_corr_masked ./ OD_corr_masked;

% ---------- Masked stats (ratio) ----------
isMasked = isnan(ratio_Flc_ODcorr);
metrics.n_masked_ratio   = sum(isMasked);
metrics.frac_masked_ratio = metrics.n_masked_ratio / numel(t);
firstValidIdx = find(~isMasked, 1, 'first');
if ~isempty(firstValidIdx)
    metrics.first_valid_ratio_time = t(firstValidIdx);
end

% ---------- 6) Lead/lag via continuous-shift NCC + derivative check ----------
validMask = ~isnan(Flc_corr_s) & ~isnan(OD_corr_s);
t_valid = t(validMask);
F_valid = Flc_corr_s(validMask);
O_valid = OD_corr_s(validMask);

if numel(t_valid) >= 5 && range(F_valid) > 0 && range(O_valid) > 0
    F0 = (F_valid - min(F_valid)); F0 = F0 / max(F0);
    O0 = (O_valid - min(O_valid)); O0 = O0 / max(O0);

    % Rising window
    Fmask = (F0 >= 0.10) & (F0 <= 0.90);
    Omask = (O0 >= 0.10) & (O0 <= 0.90);
    riseMask = Fmask | Omask;
    if nnz(riseMask) < 8
        riseMask = (F0 >= 0.05) | (O0 >= 0.05);
    end

    tW = t_valid(riseMask);
    FW = F0(riseMask);
    OW = O0(riseMask);

    good = numel(tW) >= 5 && range(FW) > 0 && range(OW) > 0;
    if good
        % Z-score in window
        FWz = (FW - mean(FW)) / std(FW);
        OWz = (OW - mean(OW)) / std(OW);

        % Continuous NCC
        Wdur = max(tW) - min(tW);
        if Wdur > 0 && numel(tW) >= 3
            dtW  = median(diff(tW));
            step = max(dtW/5, 1e-6);
            tau_grid = (-0.6*Wdur) : step : (0.6*Wdur);
            NCC = nan(size(tau_grid));
            for k = 1:numel(tau_grid)
                Oshift = interp1(tW, OWz, tW + tau_grid(k), 'linear', 'extrap');
                C = corrcoef(FWz, Oshift, 'rows', 'complete');
                NCC(k) = C(1,2);
            end
            [~,kmax] = max(NCC);
            metrics.tau_ncc = tau_grid(kmax);
        end

        % Derivative xcorr (rise-only)
        dtW = median(diff(tW));
        dFW = gradient(FW, dtW);
        dOW = gradient(OW, dtW);
        if std(dFW) > 0 && std(dOW) > 0
            dFWz = (dFW - mean(dFW)) / std(dFW);
            dOWz = (dOW - mean(dOW)) / std(dOW);
            maxLag = min(round(0.5*numel(tW)), 200);
            [xc_d,lags_d] = xcorr(dFWz, dOWz, maxLag, 'coeff');
            [~,i_d] = max(xc_d);
            metrics.tau_deriv = lags_d(i_d) * dtW;
        end
    end
end

% ----- 50% rise reference (guard division-by-zero) -----
Fspan = max(Flc_corr_s) - min(Flc_corr_s);
Ospan = max(OD_corr_s)  - min(OD_corr_s);
if Fspan > 0 && Ospan > 0
    Fnorm = (Flc_corr_s - min(Flc_corr_s)) / Fspan;
    Onorm = (OD_corr_s   - min(OD_corr_s))  / Ospan;
    t_half_Flc = interp1(Fnorm, t, 0.5, 'linear', 'extrap');
    t_half_OD  = interp1(Onorm, t, 0.5, 'linear', 'extrap');
    metrics.lead_half = t_half_OD - t_half_Flc;
end

% ---------- 13) Window for ratio slope & means ----------

%Change window over which dy/dt is calculated here



t_ratio_start = 9.5;    %%%%%%%%%% adjust as needed%%%%%%%%%%%%%%%%%
t_ratio_end   = 24;   % %%%%%%%%%%%%%adjust as needed%%%%%%%%%%%%%%%
mask_window = (t >= t_ratio_start) & (t <= t_ratio_end) & ~isnan(ratio_Flc_ODcorr);
t_fit = t(mask_window);
y_fit = ratio_Flc_ODcorr(mask_window);

% --- SG & spline means + SG peak (with bounds), guarded ---
if numel(t_fit) >= 5 && range(y_fit) > 0
    % SG parameters
    order = 2;
    frame = min(5, 2*floor(numel(t_fit)/4)+1);
    frame = max(frame, order+1 + mod(order+1,2));
    frame = min(frame, numel(t_fit) - ~mod(numel(t_fit),2));
    if mod(frame,2)==0, frame = frame-1; end
    ok_sg = frame >= (order+1) && frame <= numel(y_fit) && frame >= 3;

    if ok_sg
        y_sg = sgolayfilt(y_fit, order, frame);
        [~, g] = sgolay(order, frame);
        half = (frame-1)/2;
        dt = median(diff(t_fit));
        dy_dt = nan(size(y_fit));
        for n = (half+1):(numel(y_fit)-half)
            dy_dt(n) = dot(g(:,2), y_sg(n-half:n+half)) / dt;
        end
        metrics.avg_rate_sg_mean = mean(dy_dt, 'omitnan');
        metrics.avg_rate_sg = trapz(t_fit(~isnan(dy_dt)), dy_dt(~isnan(dy_dt))) ...
                              / (t_fit(end) - t_fit(1));

        % ---- SG peak rate + time window ----
        [peak_val, idx_peak] = max(dy_dt, [], 'omitnan');
        if ~isempty(idx_peak) && ~isnan(peak_val)
            metrics.peak_rate_sg = peak_val;
            metrics.t_peak_center = t_fit(idx_peak);
            metrics.t_peak_start  = t_fit(max(1, idx_peak - half));
            metrics.t_peak_end    = t_fit(min(numel(t_fit), idx_peak + half));
        end
    end

    % Spline mean + integral-of-derivative mean
    try
        pp  = csaps(t_fit, y_fit, 0.95);
        y0  = fnval(pp, t_fit(1));
        y1  = fnval(pp, t_fit(end));
        metrics.avg_rate_spline = (y1 - y0) / (t_fit(end) - t_fit(1));
        ppd = fnder(pp,1);
        t_dense = linspace(t_fit(1), t_fit(end), 200);
        dy_dt_spline = fnval(ppd, t_dense);
        metrics.avg_rate_spline_int = trapz(t_dense, dy_dt_spline) / (t_fit(end)-t_fit(1));
    catch
        % leave NaNs
    end
end




% ----- Return processed time-series for plotting -----
% Ensure all vectors are columns of identical length
t            = t(:);
OD_corr_s    = OD_corr_s(:);
Flc_corr_s   = Flc_corr_s(:);

metrics.t_series   = t;
metrics.OD_corr_s  = OD_corr_s;
metrics.Flc_corr_s = Flc_corr_s;

% Normalize 0–1 for shape comparison
if range(OD_corr_s) > 0
    metrics.OD_norm  = (OD_corr_s - min(OD_corr_s)) / range(OD_corr_s);
else
    metrics.OD_norm  = zeros(size(OD_corr_s));
end
if range(Flc_corr_s) > 0
    metrics.Flc_norm = (Flc_corr_s - min(Flc_corr_s)) / range(Flc_corr_s);
else
    metrics.Flc_norm = zeros(size(Flc_corr_s));
end

% 50% rise points (already computed earlier)
if exist('t_half_OD','var'),  metrics.t_half_OD  = t_half_OD;  end
if exist('t_half_Flc','var'), metrics.t_half_Flc = t_half_Flc; end



metrics = orderfields(metrics);


%%%function end 
end
