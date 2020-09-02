
filename = '../output/genBusResults.csv';
sol = csvread(filename,1,0);
ref_gen = 1;
gen_display = [1, 2, 3];
len = length(gen_display);
n = 11;
stepsize = ceil(length(sol(:,1))/20);

figure;
hold on;

for i = 1:length(gen_display)
    gen_id = gen_display(i);
    lege = ['\omega Gen ', num2str(gen_id)];
    plot(sol(:,1), sol(:,2+n*(gen_id-1)),'-','LineWidth',...
        1,'MarkerIndices',1:10:length(sol(:,1)), ...
        'DisplayName', lege);
end

hline = refline([0 1.0]);
hline.Color = 'r';
hline.LineStyle = '-.';
hline.DisplayName = 'ref. line y = 1.0';

hline = refline([0 1.005]);
hline.Color = 'blue';
hline.LineStyle = '-.';
hline.DisplayName = 'ref. line y = 1.005';

ax = gca;
ax.FontSize = 18;
ylim([0.994,1.01]);
grid on;
grid minor;
title('\omega of all generators');
xlabel('t');
ylabel('\omega values in p.u');
legend show;
legend('Location','southeast');