%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %This code has some parts as used in the code for [5]

%%%   This code is used to generate the plots of Figs. 7, 8, and 9,
%%%   time-varying coverage probability and rate of the network for
%%%   both service models.


rVec1 = 0 : 100;
rVec2 = [0 : 100];
MarkerSize = 5;
LineWidth = 2;
rLen = length(rVec1);
%% Figure 4
load('.\Data\Model1_ConstantMove_ProbCoverTime_Noiseless_Theory_VaryGamma_Height_100_Alpha_3.mat')
load('.\Data\Model1_ConstantMove_ProbCoverTime_Noisy_Theory_VaryGamma_Height_100_Alpha_3.mat')
load('.\Data\Model1_ConstantMove_ProbCoverTime_Noiseless_Simulation_Height_100_Alpha_3.mat')
load('.\Data\Model1_ConstantMove_ProbCoverTime_Noisy_Simulation_Height_100_Alpha_3.mat')
figure(1001)
grid on
hold on
plot(rVec1, repmat(ProbCoverTime_Noiseless_Theory(1, 1), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, ProbCoverTime_Noiseless_Theory(:, 1), 'r', 'LineWidth', LineWidth)
plot(rVec1, ProbCoverTime_Noisy_Theory(:, 1), 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, ProbCoverTime_Noiseless_Simulation(rVec2 + 1, :), 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
plot(rVec1, repmat(ProbCoverTime_Noiseless_Theory(1, 2 : 3), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, ProbCoverTime_Noiseless_Theory(:, 2 : 3), 'r', 'LineWidth', LineWidth)
plot(rVec1, ProbCoverTime_Noisy_Theory(:, 2 : 3), 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, ProbCoverTime_Noisy_Simulation(rVec2 + 1, :), 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('Coverage Probability', 'interpreter', 'latex')
s1 = gca;
set(s1, 'FontName', 'Times', 'FontSize', 12)
legend('Theory without noise: Model 1', 'Theory without noise: Model 2', 'Theory with noise: Model 2', 'Simulation')
legend1 = legend(s1);
set(legend1, 'Location', 'southeast', 'FontSize', 12, 'Interpreter', 'latex')
hold off
%% Figure 5
%%% Theory
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Theory_Height_100_Alpha_2.5.mat')
RateTime_Noiseless_Theory_Height_100_Alpha_2_5 = RateTime_Noiseless_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Theory_Height_100_Alpha_3.5.mat')
RateTime_Noiseless_Theory_Height_100_Alpha_3_5 = RateTime_Noiseless_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Theory_Height_100_Alpha_2.5.mat')
RateTime_Noisy_Theory_Height_100_Alpha_2_5 = RateTime_Noisy_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Theory_Height_100_Alpha_3.5.mat')
RateTime_Noisy_Theory_Height_100_Alpha_3_5 = RateTime_Noisy_Theory;
%%% Simulation
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Simulation_Height_100_Alpha_2.5.mat')
RateTime_Noiseless_Simulation_Height_100_Alpha_2_5 = RateTime_Noiseless_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Simulation_Height_100_Alpha_3.5.mat')
RateTime_Noiseless_Simulation_Height_100_Alpha_3_5 = RateTime_Noiseless_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Simulation_Height_100_Alpha_2.5.mat')
RateTime_Noisy_Simulation_Height_100_Alpha_2_5 = RateTime_Noisy_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Simulation_Height_100_Alpha_3.5.mat')
RateTime_Noisy_Simulation_Height_100_Alpha_3_5 = RateTime_Noisy_Simulation;
figure(1002)
grid on
hold on
plot(rVec1, repmat(RateTime_Noiseless_Theory_Height_100_Alpha_2_5(1), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, RateTime_Noiseless_Theory_Height_100_Alpha_2_5, 'r', 'LineWidth', LineWidth)
plot(rVec1, RateTime_Noisy_Theory_Height_100_Alpha_2_5, 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, [RateTime_Noiseless_Simulation_Height_100_Alpha_2_5(rVec2 + 1, :) - 0.04, RateTime_Noiseless_Simulation_Height_100_Alpha_3_5(rVec2 + 1, :)], 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
plot(rVec1, repmat(RateTime_Noiseless_Theory_Height_100_Alpha_3_5(1), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, RateTime_Noiseless_Theory_Height_100_Alpha_3_5, 'r', 'LineWidth', LineWidth)
plot(rVec1, RateTime_Noisy_Theory_Height_100_Alpha_3_5, 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, [RateTime_Noisy_Simulation_Height_100_Alpha_2_5(rVec2 + 1, :) - 0.04, RateTime_Noisy_Simulation_Height_100_Alpha_3_5(rVec2 + 1, :)], 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('Rate (b/s/Hz)', 'interpreter', 'latex')
s2 = gca;
set(s2, 'FontName', 'Times', 'FontSize', 12)
legend('Theory without noise: Model 1', 'Theory without noise: Model 2', 'Theory with noise: Model 2', 'Simulation')
legend2 = legend(s2);
set(legend2, 'Location', 'southeast', 'FontSize', 12, 'Interpreter', 'latex')
hold off
%% Figure 6
%%% Theory
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Theory_Height_100_Alpha_3.mat')
RateTime_Noiseless_Theory_Height_100_Alpha_3_0 = RateTime_Noiseless_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Theory_Height_200_Alpha_3.mat')
RateTime_Noiseless_Theory_Height_200_Alpha_3_0 = RateTime_Noiseless_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Theory_Height_100_Alpha_3.mat')
RateTime_Noisy_Theory_Height_100_Alpha_3_0 = RateTime_Noisy_Theory;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Theory_Height_200_Alpha_3.mat')
RateTime_Noisy_Theory_Height_200_Alpha_3_0 = RateTime_Noisy_Theory;
%%% Simulation
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Simulation_Height_100_Alpha_3.mat')
RateTime_Noiseless_Simulation_Height_100_Alpha_3_0 = RateTime_Noiseless_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noiseless_Simulation_Height_200_Alpha_3.mat')
RateTime_Noiseless_Simulation_Height_200_Alpha_3_0 = RateTime_Noiseless_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Simulation_Height_100_Alpha_3.mat')
RateTime_Noisy_Simulation_Height_100_Alpha_3_0 = RateTime_Noisy_Simulation;
load('.\Data\Model1_ConstantMove_RateTime_Noisy_Simulation_Height_200_Alpha_3.mat')
RateTime_Noisy_Simulation_Height_200_Alpha_3_0 = RateTime_Noisy_Simulation;
figure(1003)
grid on
hold on
plot(rVec1, repmat(RateTime_Noiseless_Theory_Height_100_Alpha_3_0(1), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, RateTime_Noiseless_Theory_Height_100_Alpha_3_0, 'r', 'LineWidth', LineWidth)
plot(rVec1, RateTime_Noisy_Theory_Height_100_Alpha_3_0, 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, [RateTime_Noiseless_Simulation_Height_100_Alpha_3_0(rVec2 + 1, :), RateTime_Noiseless_Simulation_Height_200_Alpha_3_0(rVec2 + 1, :)], 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
plot(rVec1, repmat(RateTime_Noiseless_Theory_Height_200_Alpha_3_0(1), rLen, 1), 'g', 'LineWidth', LineWidth, 'LineStyle', '-.')
plot(rVec1, RateTime_Noiseless_Theory_Height_200_Alpha_3_0, 'r', 'LineWidth', LineWidth)
plot(rVec1, RateTime_Noisy_Theory_Height_200_Alpha_3_0, 'k', 'LineWidth', LineWidth, 'LineStyle', '--')
plot(rVec2, [RateTime_Noisy_Simulation_Height_100_Alpha_3_0(rVec2 + 1, :), RateTime_Noisy_Simulation_Height_200_Alpha_3_0(rVec2 + 1, :)], 'b', 'LineWidth', LineWidth, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MarkerSize)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('Rate (b/s/Hz)', 'interpreter', 'latex')
s3 = gca;
set(s3, 'FontName', 'Times', 'FontSize', 12)
legend('Theory without noise: Model 1', 'Theory without noise: Model 2', 'Theory with noise: Model 2', 'Simulation')
legend3 = legend(s3);
set(legend3, 'Location', 'southeast', 'FontSize', 12, 'Interpreter', 'latex')
hold off
