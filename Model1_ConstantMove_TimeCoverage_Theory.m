 %This code has some parts as used in the code for [5]

%%%   This code is used to generate the simulation data for Figs. 7,
%%%   , time-varying coverage probability of the
%%%   network for both service models.

datetime('now')
lambda0UAV = 1e-6;
v = 45; % [km/h]
v = v / 3.6; % [m/s]
h = 100;
alpha = 3;
P = 1;
P_Edge = 0.95;
SNR_Edge = 1;
d_Edge = sqrt(-log(1 - P_Edge) / (pi * lambda0UAV) + h ^ 2);
N0 = P * d_Edge ^ (-alpha) / SNR_Edge;
tMax = 100;%600;%
dt = 1;
tVec = [0.001, dt : dt : tMax];
tLen = length(tVec);
% dg = 1;
% gamMax = 400;
% gVec = 0 : dg : gamMax;
gVec = db2pow([-6, 0, 10]);
gLen = length(gVec);
ProbCoverTime_Noiseless_Theory = zeros(tLen, gLen);
ProbCoverTime_Noisy_Theory = zeros(tLen, gLen);
it = 0;
for t = tVec
    it = it + 1;
    tic
    parfor ig = 1 : gLen
        gam = gVec(ig);
        g1 = @(ux) ux ./ (1 + 1 / gam * ((ux .^ 2 + h ^ 2) / h ^ 2) .^ (alpha / 2));
        InnerInt1 = @(u0) arrayfun(@(u0) Fun1(u0, gam, h, alpha, v, t), u0);
        InnerInt2 = @(u0) arrayfun(@(u0) Fun2(u0, gam, h, alpha, v, t), u0);
        fun01 = @(u0) exp(-2 * pi * lambda0UAV * InnerInt1(u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV));
        fun02 = @(u0) exp(-2 * pi * lambda0UAV * InnerInt2(u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV));
        fun03 = @(u0) exp(-2 * pi * lambda0UAV * InnerInt2(u0)) .* raylpdf(u0, 1 / sqrt(2 * pi * lambda0UAV)) .* exp(-gam * N0 / P * ((u0 - v * t) .^ 2 + h ^ 2) .^ (alpha / 2));
        fun1 = @(u0) arrayfun(@(u0) fun01(u0), u0);
        fun2 = @(u0) arrayfun(@(u0) fun02(u0), u0);
        fun3 = @(u0) arrayfun(@(u0) fun03(u0), u0);
        q1 = integral(fun1, 0, v * t);
        q2 = integral(fun2, v * t, inf);
        q3 = q1 * exp(-gam * N0 / P * h ^ alpha);
        q4 = integral(fun3, v * t, inf);
        q12_Exact = real(q1 + q2);
        q34_Exact = real(q3 + q4);
        ProbCoverTime_Noiseless_Theory(it, ig) = q12_Exact;
        ProbCoverTime_Noisy_Theory(it, ig) = q34_Exact;
    end
    toc
end
save(['Model1_ConstantMove_ProbCoverTime_Noiseless_Theory_VaryGamma_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'ProbCoverTime_Noiseless_Theory')
save(['Model1_ConstantMove_ProbCoverTime_Noisy_Theory_VaryGamma_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'ProbCoverTime_Noisy_Theory')
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
