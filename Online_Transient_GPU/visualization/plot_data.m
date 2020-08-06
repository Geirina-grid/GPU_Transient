clc
clear all
close all

filename = '../output/genBusResults.csv';
sol = csvread(filename,1,0);
ref_gen = 1;

gen_display = [1:3];
% gen_display = [1, 2, 3, 6, 7, 9];

% gen_display = [1:1:64];
% gen_display = 1:6;
len = length(gen_display);
mm = 3;
nn = ceil(len/mm);

n = 11;
stepsize = ceil(length(sol(:,1))/10);

%% the first figure

% hfig = figure('name', 'Generator Variable Plots');
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[0.6 mm nn mm*nn]);
% 
% for i = 1:length(gen_display)
%     subplot(mm,nn,i);
%     gen_id =  gen_display(i);
%     xlabel('t')
%     ylabel('Gen ' + string(gen_id) +' (p.u.)')
%     title('Generator ' + string(gen_id) +' plots')
%     ax = gca;
%     ax.FontSize = 12;
%     grid on;
%     grid minor;
%     
%     hold on
%     plot(sol(:,1), sol(:,8+n*(gen_id-1)),'-x','LineWidth',...
%         1,'MarkerIndices',1:stepsize:length(sol(:,1)), 'DisplayName', 'V_m');
%     plot(sol(:,1), sol(:,9+n*(gen_id-1)),'->','LineWidth',...
%         1,'MarkerIndices',1:stepsize:length(sol(:,1)), 'DisplayName', 'E_dpp');
%     plot(sol(:,1), sol(:,10+n*(gen_id-1)),'-<','LineWidth',...
%         1,'MarkerIndices',1:stepsize:length(sol(:,1)), 'DisplayName', 'E_qpp');
%     plot(sol(:,1), sol(:,11+n*(gen_id-1)),'-d','LineWidth',...
%         1,'MarkerIndices',1:stepsize:length(sol(:,1)), 'DisplayName', 'E_dp');
%     plot(sol(:,1), sol(:,12+n*(gen_id-1)),'-*','LineWidth',...
%         1,'MarkerIndices',1:stepsize:length(sol(:,1)), 'DisplayName', 'E_qp');
% end
% legend show;
% legend('Location','southeast');

%% the second figure
hfig = figure('name', 'speed and angle display');
pos = get(hfig,'position');
set(hfig,'position',pos.*[.7 1 3 3]);

subplot(1,2,1)
hold on
for i = 1:length(gen_display)
    gen_id = gen_display(i);
    lege = ['\omega Gen ', num2str(gen_id)];
    plot(sol(:,1), sol(:,2+n*(gen_id-1)),'-','LineWidth',...
        1,'MarkerIndices',1:10:length(sol(:,1)), ...
        'DisplayName', lege);
end
% hline = refline([0 1.1]);
% hline.Color = 'r';
ax = gca;
ax.FontSize = 12;
% ylim([0.99,1.01]);
grid on;
grid minor;
% hline = refline([0 1.0]);
% hline.Color = 'r';
% hline.LineStyle = '-.';
% hline.DisplayName = 'ref. line y = 1.0';
% 
% hline = refline([0 1.005]);
% hline.Color = 'blue';
% hline.LineStyle = '-.';
% hline.DisplayName = 'ref. line y = 1.005';

title('\omega of selected generators');
xlabel('t');
ylabel('\omega values in p.u');
% legend show;
% legend('Location','southeast');

% subplot(1,3,2)
% for i = 1:length(gen_display)
%     gen_id = gen_display(i);
%     hold on
%     lege = ['\delta Gen ', num2str(gen_id)];
%     plot(sol(:,1), (sol(:,5+n*(gen_id-1))) * 180/pi,'-','LineWidth',...
%         1,'MarkerIndices',1:10:length(sol(:,1)), ...
%         'DisplayName', lege);
% end
% ax = gca;
% ax.FontSize = 12;;
% grid minor;
% xlabel('t');
% title('\delta absolute of all generators');
% ylabel('\delta values in degree');
% legend show;
% legend('Location','southeast');

subplot(1,2,2)
for i = 1:length(gen_display)
    gen_id = gen_display(i);
    hold on
    lege = ['\delta Gen ', num2str(gen_id)];
    plot(sol(:,1), (sol(:,5+n*(gen_id-1)) - sol(:,5+n*(ref_gen-1))) * 180/pi,'-','LineWidth',...
        1,'MarkerIndices',1:10:length(sol(:,1)), ...
        'DisplayName', lege);
end
ax = gca;
ax.FontSize = 12;
grid on;
grid minor;
xlabel('t');
title('\delta-\delta_0 of selected generators');
ylabel('\delta values in degree');
% legend show;
% legend('Location','southeast');



% %% the third figure
% hfig = figure('name', 'Other Components Variable Plots');
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 2 3 2]);
% 
% n = 11;
% stepsize = ceil(length(sol(:,1))/20);
% 
% 
% subplot(2,2,1)
% for i = 1:length(gen_display)
%     gen_id = gen_display(i);
%     hold on
%     lege = ['\mu Gen ', num2str(gen_id)];
%     plot(sol(:,1), sol(:,3+n*(i-1)),'-','LineWidth',...
%         1,'MarkerIndices',1:10:length(sol(:,1)), ...
%         'DisplayName', lege);
% end
% ax = gca;
% ax.FontSize = 12;;
% grid on;
% grid minor;
% xlabel('t');
% title('\mu of selected generators');
% ylabel('\mu values in p.u.');
% legend show;
% legend('Location','southeast');
% 
% 
% subplot(2,2,2)
% for i = 1:len
%     gen_id = gen_display(i);
%     hold on
%     lege = ['P_T Gen ', num2str(gen_id)];
%     plot(sol(:,1), sol(:,4+n*(i-1)),'-','LineWidth', ...
%         1,'MarkerIndices',1:10:length(sol(:,1)), ...
%         'DisplayName', lege);
% end
% ax = gca;
% ax.FontSize = 12;;
% % ylim([0.99,1.01]);
% grid on;
% grid minor;
% xlabel('t');
% title('P_T of selected generators');
% ylabel('P_T values in p.u.');
% legend show;
% legend('Location','southeast');
% 
% 
% subplot(2,2,3)
% for i = 1:len
%     gen_id = gen_display(i);
%     hold on
%     lege = ['E_{fd} Gen ', num2str(gen_id)];
%     plot(sol(:,1), sol(:,6+n*(i-1)),'-','LineWidth',...
%         1,'MarkerIndices',1:10:length(sol(:,1)), ...
%         'DisplayName', lege);
% end
% ax = gca;
% ax.FontSize = 12;;
% grid on;
% grid minor;
% xlabel('t');
% title('E_{fd} of selected generators');
% ylabel('E_{fd} values in p.u.');
% legend show;
% legend('Location','southeast');
% 
% 
% subplot(2,2,4)
% for i = 1:len
%     gen_id = gen_display(i);
%     hold on
%     lege = ['V_s Gen ', num2str(gen_id)];
%     plot(sol(:,1), sol(:,7+n*(i-1)),'-','LineWidth',...
%         1,'MarkerIndices',1:10:length(sol(:,1)), ...
%         'DisplayName', lege);
% end
% ax = gca;
% ax.FontSize = 12;;
% ylim([-0.05,0.05]);
% grid on;
% grid minor;
% xlabel('t');
% title('V_s of selected generators');
% ylabel('V_s values in p.u.');
% legend show;
% legend('Location','southeast');
% 
