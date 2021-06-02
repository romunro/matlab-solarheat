%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab model for Solar Collector Group 10%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%
%%%%%%%%%%%
%Constants%
%%%%%%%%%%%
%Natural costants:
    sigma = 5.67*10^(-8);                       %Stefan-Boltzmann constant
    v_out = 4.5;                                %Velocity of air outside the solar collector
    v_in = 0;                                   %Velocity of air inside the solar collector 
%Emissivity
    e_Cu = 0.22;                                %Copper tube, value 0.052
    e_glass = 0.94;
    e_PVC = 0.91;                               %PVC tube
    e_foam = 0.9;                               %foam foil
    e_Al = 0.04;                                %Aluminum tape
    e_paint = 0.95;                             %Emissivity of black paint (for comparison)
    e_wood = 0.94;
    e_CuMatte = 0.22;                           %Emissivity of the copper tube after making it less reflective (for comparison)
%Thermal conductivity
    k_air = 0.026;                              %Air [W/(mK)]
    k_tempex =0.03;                             %Tempex plate [W/(mK)]
    k_Cu = 400;                                 %Copper tube [W/(mK)]
    k_wood = 0.12;                              %Wood casing [W/(m*K)]
    k_Po = 0.13;                                %Polyurethane tube [W/(m*K)]
    k_PVC = 0.19;                               %PVC thermal conduc [W/(mK)]
    k_foil = 0.04;                              %Polyethylene foam foil [W/(mK)]
    k_glass = 0.78;                             %Thermal conductivity of the glass plate [J/s mK]
%Heat transfer coefficient
    h_out_air = 10.45 - v_out + 10*(v_out)^(1/2);%Outside solar collector air [W/(m*K)]
    h_in_air = 10.45 - v_in + 10*(v_in)^(1/2);  %Inside solar collector air [W/(m*K)]
    h_air = 2.5;                                %Still standing 298K air [W/(m*K)]
%Specific heat
    c_water = 4186;                             %Specific heat water [J/kgK]
    c_Cu = 386;                                 %Specific heat copper [J/kgK]
    c_Al = 922;                                 %Specific heat aluminium[J/kgK]
    c_air = 1007;                               %Specific heat air[J/kgK]
%Measurements
    D_Cu = 0.012;                               %Diameter of the copper tube [m]
    D_PolyTube = 0.008;                         %Diameter of the polyurathane pump tubes [m]
    R_Cu1 = 0.005;                              %Inner radius of the copper tube [m]
    R_Cu2 = 0.006;                              %Outer radius of  the copper tube [m]
    V_air = 0.04911;                            %Volume of the air in the solar collector [m^3]
    V_Al = 8.96*10^(-5);                        %Total volume of the aluminum tape
    L_Cu_SC = 6.63;                             %Length of the copper tube in the solar collector [m]
    L_CuTube = 6.630;                           %Total length of tube that's exposed to the air [m]
    L_PolyTube1 = 3;                            %Total length of polyurathane pump tube 1 [m]
    L_PolyTube2 = 3;                            %Total length of polyurathane pump tube 2 [m]
    L_Tube_HV = 0.71;
    R_PVC2 = 0.0268;
    dx = 0.05;                                  %Thickness of the pinewood frame
%Measurements Heat storage vessel
    D_pvc = 0.050;                              %Diameter PVC tube [m]
    R_pvcThick = 0.0018;                        %PVC wall thickness
    R_pvc1 = (D_pvc/2) - R_pvcThick;            %Inner Radius PVC tube [m]
    R_pvc2 = (D_pvc/2);                         %Outer Radius PVC tube [m]
    N_insLayers = 6;                            %Number of insulating polyethylene foam foil layers [m]
    R_polyFoil = N_insLayers * 0.003;           %Radius thickness polyethylene foam foil [m]
    L_pvc = 0.71;                               %PVC length [m]
    A_HV = 2*pi*(R_pvc2+R_polyFoil)*L_pvc;      %Surface area Heat storage vessel
%Measurements for in- and outlet connectors
    D_Po = 0.012;                               %Diameter of the polyurethane tube [m]
    R_Po1 = 0.004;                              %Inner radius of the polyurethane tube [m]
    R_Po2 = 0.006;                              %Outer radius of the polyurethane tube [m]
    L_Tube_HV_to_SC = 3;                        %Half of general setup length, HV to SC [m]
    L_Tube_SC_to_HV = 3;                        %Half of general setup length, SC to HV [m]
%Areas
    A_RadCu = 4*D_Cu*1.4+3*0.5*0.25*pi*(0.0898^2-0.0778^2)+2*0.012*0.15; %Sunlit area of the copper tube
    A_RadAl = 0.670*1.4 - 4*0.012*1.4;                      %Sunlit area of the aluminum troughs
    A_AirCu = 1.4*pi*D_Cu + 3*4.9594762*10^(-3) + 0.05*pi*D_Cu + 2*0.1*pi*D_Cu; %Area of the copper tube that's exposed to the air
    A_AirTape = 1.971;                                       %Area of the aluminum tape that is exposed to the air
    A_AirSS = 2*(1.64*0.67)+2*(0.67*0.065)+2*(0.065*1.72);  %Area of solar collector that is exposed to the air
    A_frame = 2*(0.67*0.065)+2*(0.065*1.72);                %Area of sides touching pinewood
    A_vessel = 2*pi*(R_PVC2)^2 + L_Tube_HV*2*pi*R_PVC2;
%Intensity
    I_sun = 1000;                               %Intensity of the artificial sun [W/m^2]
    I_glass = 950;                              %Intensity after absorption of glass plate [W/m^2]
%Temperatures
    T_sur = 293;                                %Temperature of the surroundings [K]
    T_in = 293;                                 %Temperature of the incoming water [K]
%Density
    rho_w = 1000;                               %Density of water [kg/m^3]
    rho_Cu = 8.3*10^3;                          %Density of copper [kg/m^3]
    rho_air = 1.29;                             %Density of air [kg/m^3]
    rho_Al = 2.70*10^3;                         %Density of aluminum [kg/m^3]
%Variables
    test_length = 20;                           %Given time in test [minutes]
    flowrate = 0.1;                             %Initial flowrate [L/min]
%%
%%%%%%%%%%%%%%
%Calculations%
%%%%%%%%%%%%%%
%%Aluminium reflection and absorption rates%%
dQdt_RadAl = (1-e_Al)*A_RadAl*I_glass;                                      %Energy that is reflected by the aluminum tape [W]
dQdt_RadAl_in = (e_Al)*A_RadAl*I_glass;                                     %Energy that is absorbed by the aluminum tape [W]
dQdt_RadCu = e_Cu*A_RadCu*I_glass + e_Cu*dQdt_RadAl;                        %Energy that is absorbed by the copper tube due to radiation [W]

%%Length of simulation%%
t_end = test_length*60;                                                     %End time [s]
t_step = 0.1;                                                               %Step size [s]
Steps = t_end/t_step;                                                       %Amount of steps

%%HEAT VESSEL%%
T_HV_table = zeros(3,Steps+1);                                              %Empty vector for T_HV_table. Row 1: time t. Row 2: T_HV_table. Row 3: dQdt_SC_total.
T_HV_table(1,:) = 0:t_step:t_end;
T_HV_out=T_in;                                                              %Beginning temperature of water in the heat vessel
m_SC_water = 0.512;                                                         %Max volume of water in the solar collector
m_HV_water = 1.2;                                                           %Max volume of water in the heat vessel 
m_out_tube =1/4*pi*((D_Po -2*0.002)^2)*L_Tube_HV_to_SC ;                    %Max volume of water in the outlettubing
m_hv_frac = 0;                                                              %Initial fraction hot water to cold water

%%Thermocline HV start values%%
m_HV_new = 0;
m_HV_old = m_HV_water;
T_HV_inside = T_in;

%%Tube Heat vessel --> Solar Collector%%
m_PolyTube1_water = (1/4)*pi*D_PolyTube^2*L_PolyTube1*rho_w;                %Maximum amount of water in the tube [kg]
T_SC_in = T_in;                                                             %Starting temperature [K]

%%SOLAR COLLECTOR%%
T_SC_table = zeros(3,Steps+1);                                              %Empty vector for T_SC_table. Row 1: time t. Row 2: T_SC_table. Row 3: dQdt_SC_total.
T_SC_table(1,:) = 0:t_step:t_end;
T_SC_out=T_in;                                                              %Beginning temperature of water in the solar collector [K]

%SC glass plate calculations
d_glass = 0.004;                                                            %Thickness of the glass plate [m]
A_glass = 0.67*(1.640+0.080);                                               %Area of the glass plate [m]

%%Tube Solar Collector --> Heat vessel%%
m_PolyTube2_water = (1/4)*pi*D_PolyTube^2*L_PolyTube2*rho_w;                %Maximum amount of water in the tube [kg]
T_HV_in = T_in;                                                             %Starting temperature [K]

%%General%%
m_Cu = (L_Cu_SC*pi*R_Cu2^2 - L_Cu_SC*pi*R_Cu1^2)*rho_Cu;                    %Mass of the copper tube in the solar collector
m_air = rho_air*V_air;                                                      %Mass of the air in the solar collector
m_Al = rho_Al*V_Al;                                                         %Total mass of the aluminum tape
T_air = T_in;                                                               %Starting temperature of the air [K]
T_Al = T_in;                                                                %Starting temperature of the aluminum tape [K]
m_flow = (flowrate/60/1000)*rho_w*t_step;                                   %Mass of water inflow from pump during t_step
%%
for t=0:t_step:t_end
    %%Variable flowrate
    if (120<t && t<1080)
           xStart = 0.2;
           dx = 0.4;
           N = 8;
           pump = xStart + (0:N-1)*dx;
           flowrate = pump(round((t-60)/120));
           m_flow = (flowrate/60/1000)*rho_w*t_step;                        %Mass of water inflow from pump during t_step
    elseif (1080<t && t<1200)
           flowrate = 3;
           m_flow = (flowrate/60/1000)*rho_w*t_step;                        %Mass of water inflow from pump during t_step
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Equations for the vessel       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%Insulation heat loss%%%%
    R_HV_PVC_cond = R_pvc1/(k_PVC*(2*pi*R_pvc2*L_pvc));                     %Thermal resitance of conduction through PVC [K/W]
    R_HV_ins_Cond =  (log((R_pvc2 + R_polyFoil)/R_pvc2)) / (k_foil*A_HV);   %Thermal resitance of conduction through poly foam foil [K/W]
    R_HV_ins_Conv = 1 / (h_air*A_HV);                                       %Thermal resitance of convection outside poly foam foil [K/W]
    HL_RadOut_PVC = e_PVC * sigma * A_vessel * (T_HV_inside^4 - T_sur^4);   %Heat loss due to radiation PVC
    HL_RadOut_foam = e_foam * sigma * A_vessel * (T_HV_inside^4 - T_sur^4); %Heat loss due to radiation foam foil
    R_HV_ins_total = R_HV_PVC_cond + R_HV_ins_Cond + R_HV_ins_Conv;         %Total thermal resistance in series
    dQdt_HV_insulation = m_hv_frac*(-((T_HV_inside - T_sur ) / R_HV_ins_total)-HL_RadOut_PVC-HL_RadOut_foam);   %Total heat flow through insulation multiplied with fraction hot water
    
    %%%%%End caps heat loss%%%%
    R_HV_cap_Conv = 1 /h_air*pi*(R_pvc2+R_pvcThick).^2;                     %Thermal resistance due to convection on end caps [K/W]
    R_HV_cap1 = R_pvcThick / (k_PVC * (pi*R_pvc2).^2);                      %Thermal resistance end-cap [K/W]
    R_HV_cap2 = R_HV_cap1;                                                  %Thermal resistance same as upper end cap as the temperature in the heat storage vessel is homogeneous [K/W]
    if m_hv_frac == 1                                                       %If heatvessel is full of hot water, heat loss from both end caps is presumed
        R_HV_caps_total = R_HV_cap1 + R_HV_cap2 + 2*R_HV_cap_Conv;
        dQdt_HV_caps = -((T_HV_inside - T_sur ) / R_HV_caps_total);
    else
        R_HV_caps_total = R_HV_cap1 + R_HV_cap_Conv;                        %If heatvessel is not yet full of hot water, heat loss from only top end-cap (only 1 end-cap) is presumed
        dQdt_HV_caps = -((T_HV_inside - T_sur ) / R_HV_caps_total);
    end
    
    %%%%Total heat flow%%%%
    dQdt_HV_total = dQdt_HV_insulation + dQdt_HV_caps;                      %Total heat flow storage vessel [W/m^2]
    
    %%%%Calclate internal temperature%%%%
    delta_T = (dQdt_HV_total*t_step/(m_HV_water*c_water));                  %Resulting delta T of heat flow [K]
    T_HV_inside = T_HV_inside+delta_T;                                      %Resulting internal temperature [K]
    
    %%%%Thermocline effect%%%% 
    if (m_HV_new < m_HV_water)                                              % Cold region thermocline provides outlfow HV
        m_HV_new = m_HV_new + m_flow;
        T_HV_inside = (m_flow*T_HV_in + T_HV_inside*m_HV_new)/(m_HV_new + m_flow);
        m_hv_frac = m_HV_new / m_HV_water;                                  %Fraction of how full the HV is with 'hot water' compared to total capacity, used for fractional heat loss HV
        m_HV_old = m_HV_old - m_flow;
        T_HV_out = T_in;
    else                                                                    %Thermocline when cold water has been dispersed
        T_HV_inside = (m_flow*T_HV_in + T_HV_inside*m_HV_new)/(m_HV_new + m_flow);
        m_hv_frac = 1;                                                      %Fraction of cold water to hot water is 1, used for heat losses in HV
        T_HV_out = T_HV_inside;
    end
    
    %%%%Log HV data for plotting%%%%
    Column = round((1/t_step)*t+1);
    T_HV_table(2,Column)=T_HV_out;                                          %Assign value for T_HV_out to the right space
    T_HV_table(3,Column)=dQdt_HV_total;                                     %Assign value for dQdt_HV_total to the right space
    T_HV_table(4,Column)=T_HV_inside;                                       %Assign value for inside temperature heat vessel
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tube 1 Heat vessel --> Solar collector%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Inflow temperature%%%%
    m_PolyTube1_old = m_PolyTube1_water - m_flow;                           %Mass of water from t-t_step [kg]
    m_PolyTube1_new = m_flow;                                               %Mass of water added to the heat vessel during t_step
    T_SC_in = (m_PolyTube1_old*T_SC_in + m_PolyTube1_new*T_HV_out)/(m_PolyTube1_old+m_PolyTube1_new); %Inflow temperature [K] 
    
    %%%%Thermal resistance%%%%
    R_tube_out_Cond =  (log(R_Po2 /R_Po1)) / ((2*pi*k_Po*L_Tube_HV_to_SC ));%Thermal resitance of conduction outlet tubing [K/W]
    R_tube_out_Conv = 1 / (h_air*(2*pi*R_Po2)*L_Tube_HV_to_SC) ;           %Thermal resitance of convection outlet tubing [K/W]
    R_tube_out_total = R_tube_out_Cond + R_tube_out_Conv;                  %Total thermal resistance in series [K/W]
    
    %%%%Total heat flowflow%%%%
    dQdt_PolyTube1_total = -((T_HV_out - T_sur ) / R_tube_out_total);       %Total heat flow HV to SC[W/m^2]
    
    %%%%Outflow temperature
    delta_T = (dQdt_PolyTube1_total*t_step/(m_PolyTube1_water*c_water));    %Resulting delta T of heat flow [K]
    T_SC_in = T_SC_in+delta_T;                                              %Resulting outflow temperature [K]
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Equations for the solar collector   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_SC_old = m_SC_water - m_flow;                                         %Mass of water from t-t_step
    m_SC_new = m_flow;                                                      %Mass of water added to the solar collector during t_step
    T_SC_out = (m_SC_old*c_water*T_SC_out + m_Cu*c_Cu*T_SC_out + m_SC_new*c_water*T_SC_in)/(m_SC_old*c_water + m_Cu*c_Cu + m_SC_new*c_water);

    %%%Thermal resistances in the solar collector%%%	
    R_RadOut = e_Cu * sigma * A_AirCu * (T_SC_out^4 - T_air^4);            %Thermal resistance due to radiation [K/W]                
    R_CondOut = (log(R_Cu2/R_Cu1)) / (2*pi*k_Cu*L_CuTube);                 %Thermal resistance due to conduction in the copper tube [K/W]
    R_Cu_conv = 1 / (h_air*A_AirCu);                                       %Thermal resistance due to convection [K/W]
    R_SC_total = R_RadOut + R_CondOut + R_Cu_conv;                         %The total of all the thermal resistances [K/W]
  
    %%%%Total heat flow solar collector%%%%
    dQdt_SC_total = dQdt_RadCu - ((T_SC_out - T_air) / R_SC_total);        %Total heat flow solar collector [W/m^2]
    
    %%%%Calculate internal temperature%%%%    
    delta_T = (dQdt_SC_total*t_step)/(m_SC_water*c_water + m_Cu*c_Cu);      %Resulting delta T of heat flow [K]
    T_SC_out=T_SC_out+delta_T;                                              %Resulting outflow temperature [K]
    
    %%%%Log SC data for plotting%%%%
    Column = round((1/t_step)*t+1);
    T_SC_table(2,Column)=T_SC_out;                                          %Assign value for T to the right space
    T_SC_table(3,Column)=delta_T;                                           %Assign value for T to the right space
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tube 2 Solar collector --> Heat vessel%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_PolyTube2_old = m_PolyTube1_water - m_flow;                           %Mass of water from t-t_step
    m_PolyTube2_new = m_flow;                                               %Mass of water added to the heat vessel during t_step
   
    T_HV_in = (m_PolyTube2_old*T_HV_in + m_HV_new*T_SC_out)/(m_PolyTube2_old+m_HV_new); %Calculating the average temperature in the heat vessel, taking into account the inflow from the solar heater   
    
    R_tube_in_Cond =  (log(R_Po2 /R_Po1)) / ((2*pi*k_Po*L_Tube_SC_to_HV )); %Thermal resitance of conduction inlet tubing [K/W]
    R_tube_in_Conv = 1 / (h_air*(2*pi*R_Po2)*L_Tube_SC_to_HV) ;             %Thermal resitance of convection inlet tubing [K/W]
    R_tube_in_total = R_tube_in_Cond + R_tube_in_Conv;                      %Total thermal resistance in series [K/W]
    
    %%%%Total flow
    dQdt_PolyTube2_total = -((T_HV_in - T_sur ) / R_tube_in_total);         %Total heat flow general setup SC to HV [W/m^2]

    delta_T = (dQdt_PolyTube2_total*t_step/(m_PolyTube2_water*c_water));
    T_HV_in = T_HV_in+delta_T;
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Air Temperature           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %To calculate dQdt_air_total, first the temperature of the aluminum tape has to be determined:
    dQdt_Al_conv = h_air*A_AirTape*(T_Al-T_air);                            %Convective heat flow aluminium tape[W/m^2]
    dQdt_Al_total = dQdt_RadAl_in - dQdt_Al_conv;                           %Total heat flow aluminium tape[W/m^2]
    delta_T = (dQdt_Al_total*t_step)/(m_Al*c_Al);
    T_Al = T_Al + delta_T;                                                  %Temperature aluminium tape after delta T
    
    %Air temperature:
    dQdt_CondEnv = -(k_wood * A_frame * (T_air - T_sur))/dx;                %Heat flow to the environment through the wood by conduction
    dQdt_CondGlass = -(k_glass*A_glass*(T_air - T_sur))/d_glass;            %Heat flow to the environment through the glass by conduction
    dQdt_air_total = dQdt_Al_conv + dQdt_CondEnv + dQdt_CondGlass + dQdt_CondOut; %Total heat flow from air SC to enviroment [W/m^2]
    
    delta_T = (dQdt_air_total*t_step)/(m_air*c_air);                        %Resulting delta T of heat flow [K]
    T_air = T_air + delta_T;                                                %Resulting air temperature [K]
    
end
%%
%%%%%%%%%%
%Plotting%
%%%%%%%%%%

t_var=T_SC_table(1,:);          %Time variable
T_SC_var=T_SC_table(2,:);       %Outflow temperature solar collector
T_HV_var=T_HV_table(2,:);       %Outflow temperature heat vessel
T_HV_inside = T_HV_table(4,:);  %Inside temperature heat vessel

hold on
grid on

plot(t_var,T_SC_var);
plot(t_var,T_HV_var);
plot(t_var,T_HV_inside);
ylabel('Temperature (K)')

legend({'Outflow temperature solar collector','Outflow temperature heat vessel','Inside temperature heat vessel'}, 'Location','northwest')

xlim([0 t_end]);
xlabel('Time (s)')
title('Outflow temperatures')
saveas(gcf,'Outflow temperature.jpg')