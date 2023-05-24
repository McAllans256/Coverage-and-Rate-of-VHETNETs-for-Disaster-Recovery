%%%   This code is used to generate the simulation data for Figs
%%%   8, and 9, time-varying coverage probability  of the
%%%   network for both service models.

datetime('now')
lambda0UAV = 1e-6;
v = 45; % [km/h]
v = v / 3.6; % [m/s]
P = 1;
P_Edge = 0.95;
SNR_Edge = 1;
tMax = 300;%600;%
dt = 1;
tVec = [0.001, dt : dt : tMax];
tLen = length(tVec);
alphaVec = [2.5, 3, 3.5];
hVec = [100, 200];
for alpha = alphaVec
    for h = hVec
        if (h == 200) && ((alpha == 2.5) || (alpha == 3.5))
            continue;
        end
        disp(['alpha = ', num2str(alpha), ' and h = ', num2str(h)])
        d_Edge = sqrt(-log(1 - P_Edge) / (pi * lambda0UAV) + h ^ 2);
        N0 = P * d_Edge ^ (-alpha) / SNR_Edge;
        RateTime_Noiseless_Theory = zeros(tLen, 1);
        RateTime_Noisy_Theory = zeros(tLen, 1);
        tic
        parfor it = 1 : tLen
            t = tVec(it);
            InnerInt1 = @(gam, u0) arrayfun(@(gam, u0) Fun1(u0, gam, h, alpha, v, t), gam, u0);
            InnerInt2 = @(gam, u0) arrayfun(@(gam, u0) Fun2(u0, gam, h, alpha, v, t), gam, u0);
            fun01 = @(gam, u0) exp(-2 * pi * lambda0UAV * InnerInt1(gam, u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV)) ./ (1 + gam);
            fun02 = @(gam, u0) exp(-2 * pi * lambda0UAV * InnerInt2(gam, u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV)) ./ (1 + gam);
            fun03 = @(gam, u0) exp(-2 * pi * lambda0UAV * InnerInt1(gam, u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV)) .* exp(-gam * N0 / P * h ^ alpha) ./ (1 + gam);
            fun04 = @(gam, u0) exp(-2 * pi * lambda0UAV * InnerInt2(gam, u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV)) .* exp(-gam * N0 / P * ((u0 - v * t) .^ 2 + h ^ 2) ^ (alpha / 2)) ./ (1 + gam);
            fun1 = @(gam, u0) arrayfun(@(gam, u0) fun01(gam, u0), gam, u0);
            fun2 = @(gam, u0) arrayfun(@(gam, u0) fun02(gam, u0), gam, u0);
            fun3 = @(gam, u0) arrayfun(@(gam, u0) fun03(gam, u0), gam, u0);
            fun4 = @(gam, u0) arrayfun(@(gam, u0) fun04(gam, u0), gam, u0);
            q1 = integral2(fun1, 0, inf, 0, v * t);
            q2 = integral2(fun2, 0, inf, v * t, inf);
            q3 = integral2(fun3, 0, inf, 0, v * t);
            q4 = integral2(fun4, 0, inf, v * t, inf);
            q12_Exact = real(q1 + q2);
            q34_Exact = real(q3 + q4);
            RateTime_Noiseless_Theory(it) = q12_Exact;
            RateTime_Noisy_Theory(it) = q34_Exact;
        end
        toc
        save(['Model1_ConstantMove_RateTime_Noiseless_Theory_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'RateTime_Noiseless_Theory')
        save(['Model1_ConstantMove_RateTime_Noisy_Theory_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'RateTime_Noisy_Theory')
    end
end
% H1 = RateTime_Noiseless_Theory(1 : end - 1);
% H2 = RateTime_Noiseless_Theory(2 : end);
% H = (H1 + H2) * dt / 2;
% CumSumH = cumsum(H) ./ (dt * (1 : tLen - 1).');
% Rate_Theory = [RateTime_Noiseless_Theory(1); CumSumH];
% save('Model1_ConstantMove_Rate_Theory', 'Rate_Theory')
datetime('now')
function I = Fun1(u0, gam, h, alpha, v, t)
g1 = @(ux) ux ./ (1 + 1 / gam * ((ux .^ 2 + h ^ 2) / h ^ 2) .^ (alpha / 2));
q01 = integral(g1, 0, v * t - u0);
q02 = integral(@(ux) g1(ux) * 1 / pi .* acos((u0 ^ 2 - ux .^ 2 - v ^ 2 * t ^ 2) ./ (2 * ux * v * t)), v * t - u0, v * t + u0);
q03 = integral(g1, v * t + u0, inf);
I = q01 + q02 + q03;
end
function I = Fun2(u0, gam, h, alpha, v, t)
g2 = @(ux) ux ./ (1 + 1 / gam * ((ux .^ 2 + h ^ 2) / ((u0 - v * t) .^ 2 + h ^ 2)) .^ (alpha / 2));
q02 = integral(@(ux) g2(ux) * 1 / pi .* acos((u0 ^ 2 - ux .^ 2 - v ^ 2 * t ^ 2) ./ (2 * ux * v * t)), u0 - v * t, u0 + v * t);
q03 = integral(g2, u0 + v * t, inf);
I = q02 + q03;
end
