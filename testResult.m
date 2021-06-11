%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab plotting test results Group 10    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%%%%%%%%%%%%%
%Reading Data%
%%%%%%%%%%%%%%
%%Variables
testTableTemp = 'testingData\ogo10-temp.csv';
testTableFlow = 'testingData\ogo10-flow.csv';
numericalTableTemp = 'numericalData\numericalTemp.csv';
numericalTableFlow = 'numericalData\numericalFlow.csv';

%%Testing data
[testingTemperature] = readmatrix(testTableTemp,'VariableNamingRule','preserve');
[testingFlow] = readmatrix(testTableFlow,'VariableNamingRule','preserve');
%%Numerical data
[numericalTemperature] = readmatrix(numericalTableTemp,'VariableNamingRule','preserve');
[numericalFlow] = readmatrix(numericalTableFlow,'VariableNamingRule','preserve');

%%%%%%%%%%
%Plotting%
%%%%%%%%%%
testing_t_var = testingTemperature(:, 1)/60;          %Time variable
numerical_t_var = numericalTemperature(1, :)/60;

% %%Figure 1 for outflow temperature
figure(1);hold on
figure(1);grid on
figure(1); plot(testing_t_var, testingTemperature(:, 2)+273,'Color','blue');
figure(1); plot(numerical_t_var, numericalTemperature(2, :),':','Color','blue');
figure(1); plot(testing_t_var, testingTemperature(:, 4)+273,'Color','red');
figure(1); plot(numerical_t_var,numericalTemperature(3, :),':','Color','red');
figure(1); ylabel('Temperature (K)')
figure(1); legend({'Testing inflow temperature solar collector','Numerical inflow temperature solar collector', 'Testing outflow temperature solar collector','Numerical outflow temperature solar collector'}, 'Location','northwest')
figure(1); xlim([0 20]);
figure(1); xlabel('Time (min)')
figure(1); title('Testing result: Outflow temperatures')
figure(1); saveas(gcf,'figures\testTemperature.jpg')

%%Figure 2 for variable flowrate pump
figure(2);hold on
figure(2);grid on
figure(2); plot(testingFlow(:, 1)/3600, testingFlow(:, 2));
figure(2); plot(numericalFlow(1, :)/60, numericalFlow(2, :));
figure(2); xlim([0 20]);
figure(2); xlabel('Time (min)')
figure(2); ylabel('Flow rate (L/min)')
figure(2); title('Testing result: Variable flow rate over time')
figure(2); legend({'Testing flow rate','Numerical flow rate'}, 'Location','northwest')
figure(2); saveas(gcf,'figures\testFlow_rate.jpg')
