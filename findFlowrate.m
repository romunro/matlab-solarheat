clc;
clear all;
close all;

%% Find end temperature heat vessel based on flow rate - uses flow.m which is a copy of the main model. Flow.m needs to be updated if the main model is changed.
flowTable(1,:) = [0.1:0.1:3];
for a = flowTable
    Column = round((1/0.1)*a);      %Uses flow.m model - exactly same as Main_numerical_model.m without variable flowrate
    flowTable(2,Column) = flow(a);  %Log data for graphing
end

hold on
grid on

plot(flowTable(1,:),flowTable(2,:));
ylabel('Temperature [K]')

legend({'Heat vessel temperature'}, 'Location','southeast')
xlim([0, 3]);
xlabel('Flow rate [L/min]')
title('Heat vessel end temperatures depending on flow rate')
saveas(gcf,'findFlowrate.jpg')