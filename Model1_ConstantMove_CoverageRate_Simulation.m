                    %%%
                    %This code has some parts as used in the code for [5]
%%%                                                                     %%%
%%%   This code is used to generate the simulation data for Figures. 7,    %%%
%%%   8, and 9, time-varying coverage probability and rate of the       %%%
%%%   network for both service models.                                  %%%


pc = parcluster('local');
parpool(pc, str2double(getenv('SLURM_CPUS_ON_NODE')));

datetime('now')
lambda0UAV = 1e-6;
R_UAV = 1e5;
NumUAV_Initial = lambda0UAV * pi * R_UAV ^ 2;
v = 45; % [km/h]
v = v / 3.6; % [m/s]
P = 1;
P_Edge = 0.95;
SNR_Edge = 1;
dt = 1;%10;
tMax = 300;%600;%
tVec = dt : dt : tMax;
tLen = length(tVec);
% dg = 1;
% gamMax = 400;
% gVec = 0 : dg : gamMax;
gVec = [db2pow(-6), 1, 10];
gLen = length(gVec);
Realizations = 1e5;%1e6;
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
        ProbCoverTime_Noiseless_Simulation = zeros(tLen + 1, gLen);
        ProbCoverTime_Noisy_Simulation = zeros(tLen + 1, gLen);
        SINR_Noiseless = zeros(Realizations, tLen + 1);
        SINR_Noisy = zeros(Realizations, tLen + 1);
        tic
        parfor i = 1 : Realizations
            NumUAV = poissrnd(NumUAV_Initial);
            PosUAV_Range = unifrnd(0, 1, NumUAV, 1);
            PosUAV_Range = R_UAV * sqrt(PosUAV_Range);
            PosUAV_Theta = unifrnd(0, 2 * pi, NumUAV, 1);
            PosUAV = repmat(PosUAV_Range, 1, 2) .* [cos(PosUAV_Theta), sin(PosUAV_Theta)];
            [~, u0Ind] = min(PosUAV_Range);
            RecoverTheta = pi + PosUAV_Theta(u0Ind);
            Temp_Noiseless = zeros(1, tLen + 1);
            Temp_Noisy = zeros(1, tLen + 1);
            DisplacedPosUAV = PosUAV;
            DisplacedTheta = unifrnd(0, 2 * pi, NumUAV, 1);
            DisplacedTheta(u0Ind) = RecoverTheta;
            vd = v * dt * [cos(DisplacedTheta), sin(DisplacedTheta)];
            PosUAV3D = [DisplacedPosUAV, h * ones(NumUAV, 1)];
            DistanceUAV2UE = sqrt(sum(PosUAV3D .^ 2, 2));
            FadingUAV2UE = exprnd(1, NumUAV, tLen + 1);
            TotalPower = sum(P * FadingUAV2UE(:, 1) .* (DistanceUAV2UE .^ (-alpha)));
            ReceivedPower = P * FadingUAV2UE(u0Ind, 1) * (DistanceUAV2UE(u0Ind) ^ (-alpha));
            Interference = TotalPower - ReceivedPower;
            Temp_Noiseless(1) = ReceivedPower / Interference;
            Temp_Noisy(1) = ReceivedPower / (N0 + Interference);
            it = 1;
            for t = tVec
                it = it + 1;
                if norm(DisplacedPosUAV(u0Ind, :)) < v * dt
                    DisplacedPosUAV = DisplacedPosUAV + vd;
                    DisplacedPosUAV(u0Ind, :) = [0, 0];
                else
                    DisplacedPosUAV = DisplacedPosUAV + vd;
                end
                PosUAV3D = [DisplacedPosUAV, h * ones(NumUAV, 1)];
                DistanceUAV2UE = sqrt(sum(PosUAV3D .^ 2, 2));
                % FadingUAV2UE = exprnd(1, NumUAV, 1);
                TotalPower = P * FadingUAV2UE(:, it).' * (DistanceUAV2UE .^ (-alpha));
                ReceivedPower = P * FadingUAV2UE(u0Ind, it) * (DistanceUAV2UE(u0Ind) ^ (-alpha));
                Interference = TotalPower - ReceivedPower;
                Temp_Noiseless(it) = ReceivedPower / Interference;
                Temp_Noisy(it) = ReceivedPower / (N0 + Interference);
            end
            SINR_Noiseless(i, :) = Temp_Noiseless;
            SINR_Noisy(i, :) = Temp_Noisy;
        end
        toc
        for it = 1 : tLen + 1
            for ig = 1 : gLen
                gam = gVec(ig);
                ProbCoverTime_Noiseless_Simulation(it, ig) = nnz(SINR_Noiseless(:, it) >= gam) / Realizations;
                ProbCoverTime_Noisy_Simulation(it, ig) = nnz(SINR_Noisy(:, it) >= gam) / Realizations;
            end
        end
        RateTime_Noiseless_Simulation = mean(log(1 + SINR_Noiseless)).';
        RateTime_Noisy_Simulation = mean(log(1 + SINR_Noisy)).';
        save(['Model1_ConstantMove_ProbCoverTime_Noiseless_Simulation_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'ProbCoverTime_Noiseless_Simulation')
        save(['Model1_ConstantMove_ProbCoverTime_Noisy_Simulation_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'ProbCoverTime_Noisy_Simulation')
        save(['Model1_ConstantMove_RateTime_Noiseless_Simulation_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'RateTime_Noiseless_Simulation')
        save(['Model1_ConstantMove_RateTime_Noisy_Simulation_Height_', num2str(h), '_Alpha_', num2str(alpha)], 'RateTime_Noisy_Simulation')
    end
end
datetime('now')
