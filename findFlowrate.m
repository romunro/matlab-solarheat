clc;
clear all;

flowTable(1,:) = [0.1:0.1:3];

for a = flowTable
    Column = round((1/0.1)*a);
    flowTable(2,Column) = flow(a);
end

hold on
grid on

plot(flowTable(1,:),flowTable(2,:));
ylabel('Temperature [K]')

legend({'Heat vessel temperature'}, 'Location','southeast')

xlim([0, 3]);
xlabel('Flow rate [L/min]')
title('Heat vessel end temperatures depending on flow rate')