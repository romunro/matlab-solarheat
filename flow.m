function y = flow(flowrate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab model for Solar Collector %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%
%Constants%
%%%%%%%%%%%

%Natural costants:
    sigma = 5.67*10^(-8);                       %Stefan-Boltzmann constant
%Convection
    h_air = 2.5;                                    %Convective heat transfer Coefficient h. Air, free 2.5 - 25 W/(m^2*k)
%Emissivity
    e_Cu = 0.22;                                %Copper tube, value 0.052
    e_Al = 0.04;                                %Aluminum tape
    e_paint = 0.95;                             %Emissivity of black paint (for comparison)
    e_CuMatte = 0.22;                           %Emissivity of the copper tube after making it less reflective (for comparison)
%Thermal conductivity
    k_air = 0.026;                              %Air [W/(mK)]
    k_tempex =0.03;                             %Tempex plate [W/(mK)]
    k_Cu = 400;                                 %Copper tube [W/(mK)]
    k_Po = 0.13;                                %Polyurethane tube [W/(m*K)]
    k_pvc = 0.19;                               %PVC thermal conduc [W/(mK)]
    k_polyFoil = 0.04;                          %Polyethylene foam foil [W/(mK)]
%Specific heat
    c_water = 4186;                             %[J/kgK]
    c_Cu = 386;                                 %[J/kgK]
%Measurements
    D_Cu = 0.012;                               %Diameter of the copper tube [m]
    D_PolyTube = 0.008;                         %Diameter of the polyurathane pump tubes [m]
    R_Cu1 = 0.005;                              %Inner diameter of the copper tube [m]
    R_Cu2 = 0.006;                              %Outer diameter of  the copper tube [m]
    L_TubeAir = 5.81973313;                     %Total length of tube that's exposed to the air [m]
    L_TubeTempex = 0.4446626;                   %Total length of tube that's exposed to the tempex [m]
    L_PolyTube1 = 3;                            %Total length of polyurathane pump tube 1 [m]
    L_PolyTube2 = 3;                            %Total length of polyurathane pump tube 2 [m]
%Measurements Heat storage vessel
    D_pvc = 0.050;                              %Diameter PVC tube [m]
    R_pvcThick = 0.0018;                        %PVC wall thickness
    R_pvc1 = (D_pvc/2) - R_pvcThick;                %Inner Radius PVC tube [m]
    R_pvc2 = (D_pvc/2);                         %Outer Radius PVC tube [m]
    N_insLayers = 6;                            %Number of insulating polyethylene foam foil layers [m]
    R_polyFoil = N_insLayers * 0.003;           %Radius thickness polyethylene foam foil [m]
    L_pvc = 0.71;                               %PVC length [m]
    A_HV = 2*pi*(R_pvc2+R_polyFoil)*L_pvc;      %Surface area Heat storage vessel
%Measurements for in- and outlet connectors
    D_Po = 0.012;                               %Diameter of the polyurethane tube [m]
    R_Po1 = 0.004;                              %Inner radius of the polyurethane tube [m]
    R_Po2 = 0.006;                              %Outer radius of the polyurethane tube [m]
    L_Tube_HV_to_SC = 3;                        %Estimated length of tube that goes from the heat storage vessel to the solar collector [m]
    L_Tube_SC_to_HV = 3;                        %Estimated length of tube that goes from the solar collector to the heat storage vessel[m]]
%Areas
    A_RadCu = 4*D_Cu*1.4+3*0.5*0.25*pi*(0.0898^2-0.0778^2)+2*0.012*0.15; %Sunlit area of the copper tube
    A_RadAl = 0.670*1.4 - 4*0.012*1.4;          %Sunlit area of the aluminum troughs
    A_AirCu = 1.4*pi*D_Cu + 3*4.9594762*10^(-3) + 0.05*pi*D_Cu + 2*0.1*pi*D_Cu; %Area of the copper tube that's exposed to the air
%Intensity
    I_sun = 1000;                               %Intensity of the artificial sun [W/m^2]
%Temperatures
    T_sur = 293;                                %Temperature of the surroundings [K]
    T_in = 293;                                 %Temperature of the incoming water [K]
%Density
    rho_w = 1000;                               %Density of water [kg/m^3]
%Variables
    t_end = 20*60;                              %End time [s]
    t_step = 0.1;                               %Step size, do not change [s]
    m_SC_water = 0.512;                         %Max mass of water in the solar collector [kg]
    m_HV_water = 1.2;                           %Max mass of water in the heat vessel [kg]
    %flowrate = massFlow;                             %Value between 0.1-3.0 [L/min]

%%

%%%%%%%%%%%%%%
%Calculations%
%%%%%%%%%%%%%%

%Intensity after passing through the glass:
I_glass = 970; %Don't know yet, can't figure it out. Internet says around 2-4% is absorbed

dQdt_RadAl = (1-e_Al)*A_RadAl*I_glass;                                      %Energy that is reflected by the aluminum tape [W]
dQdt_RadCu = e_Cu*A_RadCu*I_glass + e_Cu*dQdt_RadAl;                        %Energy that is absorbed by the copper tube due to radiation [W]

Steps = t_end/t_step;                                                       %Amount of steps

%%HEAT VESSEL%%
T_HV_table = zeros(3,Steps+1);                                              %Empty vector for T_HV_table. Row 1: time t. Row 2: T_HV_table. Row 3: dQdt_SC_total.
T_HV_table(1,:) = 0:t_step:t_end;
T_HV_out=T_in;                                                              %Beginning temperature of water in the heat vessel
m_SC_water = 0.512; %Max volume of water in the solar collector
m_HV_water = 1.2; %Max volume of water in the heat vessel %%%CHECK VALUE%%%%
m_in_tube =1/4*pi*((D_Po -2*0.002)^2)*L_Tube_SC_to_HV ; %Max volume of water in the inlet tubing %%%CHECK VALUE%%%%
m_out_tube =1/4*pi*((D_Po -2*0.002)^2)*L_Tube_HV_to_SC ; %Max volume of water in the outlettubing %%%CHECK VALUE%%%%

%%Tube Heat vessel --> Solar Collector%%
m_PolyTube1_water = (1/4)*pi*D_PolyTube^2*L_PolyTube1*rho_w;                %Maximum amount of water in the tube [kg]
T_SC_in = T_in;                                                             %Starting temperature [K]

%%SOLAR COLLECTOR%%
T_SC_table = zeros(3,Steps+1);                                              %Empty vector for T_SC_table. Row 1: time t. Row 2: T_SC_table. Row 3: dQdt_SC_total.
T_SC_table(1,:) = 0:t_step:t_end;
T_SC_out=T_in;                                                              %Beginning temperature of water in the solar collector

%%Tube Solar Collector --> Heat vessel%%
m_PolyTube2_water = (1/4)*pi*D_PolyTube^2*L_PolyTube2*rho_w;                %Maximum amount of water in the tube [kg]
T_HV_in = T_in;                                                             %Starting temperature [K]

%%General%%
m_flow = (flowrate/60/1000)*rho_w*t_step;                                   %Mass of water inflow from pump during t_step

%%
for t=0:t_step:t_end
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Equations for the vessel       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_HV_old = m_HV_water - m_flow;                                         %Mass of water from t-t_step
    m_HV_new = m_flow;                                                      %Mass of water added to the heat vessel during t_step
    T_HV_out = (m_HV_old*T_HV_out + m_HV_new*T_HV_in)/(m_HV_old+m_HV_new);  %Calculating the average temperature in the heat vessel, taking into account the inflow from the solar heater
    
    %%%%%Insulation heat loss 
    R_HV_ins_PVC = R_pvc1/(k_pvc*(2*pi*R_pvc2*L_pvc));
    R_HV_ins_Cond =  (log((R_pvc2 + R_polyFoil)/R_pvc2)) / (k_polyFoil*A_HV);   %Thermal resitance of conduction through poly foam foil [K/W]
    R_HV_ins_Conv = 1 / (h_air*A_HV);                                           %Thermal resitance of convection outside poly foam foil [K/W]
    R_HV_ins_total = R_HV_ins_PVC + R_HV_ins_Cond + R_HV_ins_Conv;              %Total thermal resistance in series
    dQdt_HV_insulation = -((T_HV_out - T_sur ) / R_HV_ins_total);               %Total heat flow through insulation
    
    %%%%%End caps heat loss
    R_HV_cap_Conv = 1 /h_air*pi*(R_pvc2+R_pvcThick).^2;                         %Thermal resistance due to convection on end caps
    R_HV_cap1 = R_pvcThick / (k_pvc * (pi*R_pvc2).^2);                          
    R_HV_cap2 = R_HV_cap1;                                                      %Thermal resistance same as upper end cap as the temperature in the heat storage vessel is homogeneous
    R_HV_caps_total = R_HV_cap1 + R_HV_cap2 + 2*R_HV_cap_Conv;
    dQdt_HV_caps = -((T_HV_out - T_sur ) / R_HV_caps_total);
    
    
    %%%%Total heatflow
    dQdt_HV_total = dQdt_HV_insulation + dQdt_HV_caps;                          %Total heat flow storage vessel
    
    delta_T = (dQdt_HV_total*t_step/(m_HV_water*c_water));
    T_HV_out = T_HV_out+delta_T;
    
    Column = round((1/t_step)*t+1);
    T_HV_table(2,Column)=T_HV_out;                                          %Assign value for T_HV_out to the right space
    T_HV_table(3,Column)=dQdt_HV_total;                                     %Assign value for dQdt_HV_total to the right space
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tube 1 Heat vessel --> Solar collector%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_PolyTube1_old = m_PolyTube1_water - m_flow;                           %Mass of water from t-t_step
    m_PolyTube1_new = m_flow;                                               %Mass of water added to the heat vessel during t_step
    T_SC_in = (m_PolyTube1_old*T_SC_in + m_HV_new*T_HV_out)/(m_PolyTube1_old+m_PolyTube1_new); %Calculating the average temperature in the heat vessel, taking into account the inflow from the solar heater   
    
    R_tube_out_Cond =  (log(R_Po2 /R_Po1)) / ((2*pi*k_Po*L_Tube_HV_to_SC )); %Thermal resitance of conduction outlet tubing [K/W]
    R_tube_out_Conv = 1 / (h_air*(2*pi*R_Po2)*L_Tube_HV_to_SC) ;             %Thermal resitance of convection outlet tubing [K/W]
    R_tube_out_total = R_tube_out_Cond + R_tube_out_Conv;                      %Total thermal resistance in series
    
    dQdt_PolyTube1_total = -((T_HV_in - T_sur ) / R_tube_out_total);                    %Heat loss of conduction outlet tubing [W]
    
    delta_T = (dQdt_PolyTube1_total*t_step/(m_PolyTube1_water*c_water));
    T_SC_in = T_SC_in+delta_T;
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Equations for the solar collector   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_SC_old = m_SC_water - m_flow;                                         %Mass of water from t-t_step
    m_SC_new = m_flow;                                                      %Mass of water added to the solar collector during t_step
    T_SC_out = (m_SC_old*T_SC_out + m_SC_new*T_SC_in)/(m_SC_old+m_SC_new);  %Calculating the average temperature in the solar collector, taking into account the inflow from the vessel
    
    dQdt_RadOut = -e_Cu*sigma*A_AirCu*(T_SC_out^4 - T_sur^4);               %Heat loss due to radiation (assuming that T_sur remains constant)
    dQdt_CondOut = 0;                                                       %Heat loss due to conduction
    dQdt_ConvOut = 0;                                                       %Heat loss due to convection
    dQdt_SC_total = dQdt_RadCu + dQdt_RadOut + dQdt_CondOut;                %Total energy
    
    delta_T = (dQdt_SC_total*t_step/(m_SC_water*c_water));
    T_SC_out=T_SC_out+delta_T;
    
    Column = round((1/t_step)*t+1);
    T_SC_table(2,Column)=T_SC_out;                                          %Assign value for T to the right space
    T_SC_table(3,Column)=dQdt_SC_total;                                     %Assign value for T to the right space
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tube 2 Solar collector --> Heat vessel%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m_PolyTube2_old = m_PolyTube1_water - m_flow;                                                   %Mass of water from t-t_step
    m_PolyTube2_new = m_flow;                                                                       %Mass of water added to the heat vessel during t_step
    T_HV_in = (m_PolyTube2_old*T_HV_in + m_HV_new*T_SC_out)/(m_PolyTube2_old+m_PolyTube1_new);      %Calculating the average temperature in the heat vessel, taking into account the inflow from the solar heater
    
    R_tube_in_Cond =  (log(R_Po2 /R_Po1)) / ((2*pi*k_Po*L_Tube_SC_to_HV )); %Thermal resitance of conduction inlet tubing [K/W]
    R_tube_in_Conv = 1 / (h_air*(2*pi*R_Po2)*L_Tube_SC_to_HV) ;             %Thermal resitance of convection inlet tubing [K/W]
    R_tube_in_total = R_tube_in_Cond + R_tube_in_Conv;                      %Total thermal resistance in series
    
    dQdt_PolyTube2_total = -((T_HV_in - T_sur ) / R_tube_in_total);                    %Heat loss of conduction inlet tubing [W]
    
    delta_T = (dQdt_PolyTube2_total*t_step/(m_PolyTube2_water*c_water));
    T_HV_in = T_HV_in+delta_T;
end
y = max(T_HV_table(2,:));
end