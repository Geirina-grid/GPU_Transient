close all
clc 
clear all

filename = '../output/allBusResults.csv';
buses = csvread(filename,1,0);
numBus = 9;
step = floor(numBus/1);
load_display = 1:9;
% load_display = [1:step:numBus, 200];

% load_display = [4, 5, 7, 9, 10, 11, 12, 13, 14];
% load_display = [1, 10, 20, 30, 40, 50, 70, 80, 100, 120, 150, 170, 172];

hold on
for i = 1:length(load_display)
    id = load_display(i)+1;
    lege = ['V_m Bus ', num2str(id-1)];
    plot(buses(:,1), sqrt(buses(:,id).^2 + buses(:, numBus+id).^2), ...
        '-','LineWidth', 1, 'MarkerIndices', 1:10:length(buses(:,1)), ...
        'DisplayName', lege);
    hold on;
end

% xlim([9.8 14]);

% hline = refline([0 1.06]);
% hline.Color = 'r';
% hline.LineStyle = '-.';
% hline.DisplayName = 'ref. line y = 1.06';
% 
% hline = refline([0 1.045]);
% hline.Color = 'r';
% hline.LineStyle = '-.';
% hline.DisplayName = 'ref. line y = 1.045';

ax = gca;
ax.FontSize = 12;
grid on;
grid minor;

% importdata("d_009_fault_mdl_01.m")
% hold on;
% for i = 2:4
%     var_legend(i-1)
%     plot(output_data(:,1), output_data(:,i), '-.','LineWidth',1, ...
%         'DisplayName', char(var_legend(i-1)));
% end
% ylim([1.025, 1.065]);
% 

xlabel('t');
ylabel('Voltage (p.u.)');
legend show;
legend('Location','southeast');

