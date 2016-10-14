close all; clear all;

SimTime = 500;
num_cstate = 4;
num_q = 3;
deltaT = 0.2;

MatlabSimfile_APMPC_Discrete = fopen('OutFiles/MatlabSimfile_Discrete.txt', 'r');
SimData_apmpc_discrete = fscanf(MatlabSimfile_APMPC_Discrete, '%f', [9 inf]);
fclose(MatlabSimfile_APMPC_Discrete);
q_apmpc_sim = SimData_apmpc_discrete(1, :)';
x_apmpc_sim = SimData_apmpc_discrete(2:1+num_cstate, :)';
b_apmpc_sim = SimData_apmpc_discrete((num_cstate+2):(1+num_cstate+num_q), :);
control_apmpc = SimData_apmpc_discrete(end, :)';
time_apmpc = 0:deltaT:(deltaT*(length(q_apmpc_sim)-1));


MatlabSimfile_MPC_Discrete = fopen('./MPCComparison/OutFiles/MatlabSimfile_MPC_Discrete.txt', 'r');
SimData_mpc_discrete = fscanf(MatlabSimfile_MPC_Discrete, '%f', [6 inf]);
fclose(MatlabSimfile_MPC_Discrete);
q_mpc_sim = SimData_mpc_discrete(1, :)';
x_mpc_sim = SimData_mpc_discrete(2:5, :)';
control_mpc = SimData_mpc_discrete(6, :)';
time_mpc = 0:deltaT:(deltaT*(length(q_mpc_sim)-1));

mode = {'Oblivious', 'Aggressive', 'Courteous', 'Reasonable'};
driver_apmpc_mode = mode{q_apmpc_sim(1)}
driver_mpc_mode = mode{floor(q_mpc_sim(1)/3)+1}


%%
figureposition = [0 500 2048 300]; figurelinewidth = 1.5; figurefontsize = 40; markersize = 8;
savedir = './figures/';

% draw (x_t(1) - x_t(3))
figure('Position', figureposition)
plot(time_apmpc, x_apmpc_sim(:, 1) - x_apmpc_sim(:, 3), 'k', 'LineWidth', figurelinewidth+0.5); hold on;
plot(time_mpc, x_mpc_sim(:, 1) - x_mpc_sim(:, 3), 'b', 'LineWidth', figurelinewidth+0.5);
%plot(x_sim(:, 3), 'b', 'LineWidth', figurelinewidth+0.5); hold on;
xlabel('t', 'FontSize', figurefontsize);
ylabel('|d_h - d_r|', 'FontSize', figurefontsize);
%title('x(t)', 'FontSize', figurefontsize);hold on;
legend('Our method', 'Normal MPC', 'Location', [.85,.7,.1,.2]);
set(gca, 'FontSize', figurefontsize);
savepath = strcat(savedir, 'positions_APMPC_MPC_', driver_apmpc_mode, '.png');
set(gcf,'PaperPositionMode','auto'); print(gcf, '-dpng', savepath);
