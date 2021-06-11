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
% testingTemperature.Properties.VariableNames
%%%%%%%%%%
%Plotting%
%%%%%%%%%%
testing_t_var = testingTemperature(:, 1)/60;          %Time variable
numerical_t_var = numericalTemperature(1, :)/60;
% T_SC_var=testingTemperature(:, 2);       %Outflow temperature solar collector

% T_HV_var=T_HV_table(2,:);       %Outflow temperature heat vessel
% T_HV_inside = T_HV_table(4,:);  %Inside temperature heat vessel
% t_var_flow = Flowrate_Table(1,:)/60;
% V_flowrate = Flowrate_Table(2,:); %Variable flowrate for pump
% 
% 
% %%Figure 1 for outflow temperature
figure(1);hold on
figure(1);grid on
figure(1); plot(testing_t_var, testingTemperature(:, 2)+273,'Color','blue');
figure(1); plot(numerical_t_var, numericalTemperature(2, :),':','Color','blue');
figure(1); plot(testing_t_var, testingTemperature(:, 4)+273,'Color','red');
figure(1); plot(numerical_t_var,numericalTemperature(3, :),':','Color','red');
% 
% figure(1); annotation('textarrow',[0.8 0.9], [0.82 0.78] ,'String','T = 314.6  K ');
% figure(1); annotation('textarrow', [0.45 0.33], [0.25 0.25], 'String', 'Thermocline effect');
% 
figure(1); ylabel('Temperature (K)')
figure(1); legend({'Testing inflow temperature solar collector','Numerical inflow temperature solar collector', 'Testing outflow temperature solar collector','Numerical outflow temperature solar collector'}, 'Location','northwest')

figure(1); xlim([0 20]);
figure(1); xlabel('Time (min)')
figure(1); title('Outflow temperatures')
figure(1); saveas(gcf,'figures\testTemperature.jpg')

%%Figure 2 for variable flowrate pump
% figure(2);hold on
% figure(2);grid on
% figure(2); plot(t_var_flow,V_flowrate);
% figure(2); xlim([0 t_end/60]);
% figure(2); xlabel('Time (min)')
% figure(2); ylabel('Flow rate (L/min)')
% figure(2); title('Variable flow rate over time')
% figure(2); legend({'Flow rate pump'}, 'Location','northwest')
% figure(2); saveas(gcf,'Flow_rate.jpg')
