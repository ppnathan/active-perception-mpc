close all; clear all;

SimTime = 500;
num_cstate = 4;
num_q = 3;

MatlabSimfile_Discrete = fopen('OutFiles/MatlabSimfile_Discrete.txt', 'r');
SimData_discrete = fscanf(MatlabSimfile_Discrete, '%f', [9 inf]);
fclose(MatlabSimfile_Discrete);
q_sim = SimData_discrete(1, :)';
x_sim = SimData_discrete(2:1+num_cstate, :)';
b_sim = SimData_discrete((num_cstate+2):(1+num_cstate+num_q), :);
control = SimData_discrete(end, :)';

mode = {'oblivious', 'aggressive', 'courteous'};%, 'reasonable'};
driver_mode = q_sim(1);


%%
figureposition = [0 500 2048 300]; figurelinewidth = 1.5; figurefontsize = 40; markersize = 8;
savedir = './figures/';

% draw (x_t(1) - x_t(3))
figure('Position', figureposition)
plot(x_sim(:, 1) - x_sim(:, 3), 'k', 'LineWidth', figurelinewidth+0.5); hold on;
%plot(x_sim(:, 3), 'b', 'LineWidth', figurelinewidth+0.5); hold on;
xlabel('t', 'FontSize', figurefontsize);
ylabel('Position', 'FontSize', figurefontsize);
%title('x(t)', 'FontSize', figurefontsize);hold on;
%legend('x_1(t)', 'x_2(t)', 'Location', [.9,.7,.1,.2]);
set(gca, 'FontSize', figurefontsize);
savepath = [savedir 'positions_MPC_' driver_mode '.png'];
set(gcf,'PaperPositionMode','auto'); print(gcf, '-dpng', savepath);

% draw (x_t(2), x_t(4))
figure('Position', figureposition)
plot(x_sim(:, 2), 'b', 'LineWidth', figurelinewidth+0.5); hold on;
plot(x_sim(:, 4), 'k', 'LineWidth', figurelinewidth+0.5);
xlabel('t', 'FontSize', figurefontsize); ylabel('Velocities', 'FontSize', figurefontsize);
%title('v(t)', 'FontSize', figurefontsize);hold on;
legend('v_1(t)', 'v_2(t)', 'Location', [.9,.7,.1,.2]);
set(gca, 'FontSize', figurefontsize);
savepath = [savedir 'velocities_MPC_' driver_mode '.png'];
set(gcf,'PaperPositionMode','auto'); print(gcf, '-dpng', savepath);

%draw P(t);
figure('Position', figureposition)
plot(b_sim(1, :), 'k', 'LineWidth', figurelinewidth); hold on;
plot(b_sim(2, :), 'b', 'LineWidth', figurelinewidth);
plot(b_sim(3, :), 'r', 'LineWidth', figurelinewidth);
xlabel('t', 'FontSize', figurefontsize); ylabel('P_t', 'FontSize', figurefontsize-5);
%title('P(q=0, t)', 'FontSize', figurefontsize);hold on;
legend('Oblivious', 'Aggressive', 'Courteous', 'Location', [.85,.6,.1,.2]);
set(gca,'YLim',[0 1], 'FontSize', figurefontsize);
savepath = [savedir 'P_t_HMSHS_' driver_mode '.png'];
set(gcf,'PaperPositionMode','auto'); print(gcf, '-dpng', savepath);

% draw sigma(t)
figure('Position', figureposition)
plot(control, 'b', 'LineWidth', figurelinewidth);
xlabel('t', 'FontSize', figurefontsize); ylabel('sigma_1', 'FontSize', figurefontsize);
%title('\sigma (t)', 'FontSize', figurefontsize);hold on;
% set(gca,'YLim',[-0.5 2.5], 'FontSize', figurefontsize, ...
%     'YTick',[0 1 2], 'YTickLabel',['u = 0 '; 'u = 1 '; 'u = -1']);
savepath = [savedir 'sigmas_MPC_' driver_mode '.png'];
set(gcf,'PaperPositionMode','auto'); print(gcf, '-dpng', savepath);


%%
if(1)
    
figure('Position', figureposition)
change_lane_pos = -50;
left_most_pos = -510;
savepath = [savedir 'LaneMergingSim_MPC_' driver_mode '.avi'];
LaneMergingVideo = VideoWriter(savepath);
LaneMergingVideo.FrameRate = 20;
open(LaneMergingVideo);
for i = 1:size(x_sim, 1)
    plot([left_most_pos, change_lane_pos], [0, 0], 'k', 'LineWidth', figurelinewidth + 3); hold on;
    plot([left_most_pos, 50], [2, 2], 'k', 'LineWidth', figurelinewidth + 3);
    plot([left_most_pos, change_lane_pos], [-2, -2], 'k', 'LineWidth', figurelinewidth + 3);
    plot([change_lane_pos, 0], [-2, 0], 'k', 'LineWidth', figurelinewidth + 3);
    plot([0, 50], [0, 0], 'k', 'LineWidth', figurelinewidth + 3);
    title(['A ' mode{floor(q_sim(1)/3)+1} ' driver vs MPC.'], 'FontSize', figurefontsize);
    set(gca,'XLim', [left_most_pos 50], 'FontSize', figurefontsize);
    length_car = 7;
    width_car = 0.4;
    
    rectangle('Position', [x_sim(i, 1) - length_car*0.5, 1 - width_car*0.5, length_car, width_car], ...
        'FaceColor', 'b');
    if(x_sim(i, 3) < change_lane_pos)
        rectangle('Position', [x_sim(i, 3) - length_car*0.5, -1 - width_car*0.5, length_car, width_car], ...
            'FaceColor', 'k');
    elseif(x_sim(i, 3) >= change_lane_pos && x_sim(i, 3) < 0)
        y_tmp = - 2 / change_lane_pos * (x_sim(i, 3) - change_lane_pos/2);
        rectangle('Position', [x_sim(i, 3) - length_car*0.5, y_tmp - width_car*0.5, length_car, width_car], ...
            'FaceColor', 'k');
    else 
        rectangle('Position', [x_sim(i, 3) - length_car*0.5, 1 - width_car*0.5, length_car, width_car], ...
            'FaceColor', 'k');
    end
    hold off;
    frame = getframe(gcf);
    writeVideo(LaneMergingVideo, frame);
%     pause(1/LaneMergingVideo.FrameRate);
end

close(LaneMergingVideo);
 
    
end