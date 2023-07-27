classdef SimulationParameters
    properties
        %no default

        %depth [m], positive
        depth = 6000
        %mass flowrate [kg/s], positive
        m_dot;
        %values with default
        %inside well radius [m] ~18 3/8"
        well_radius = 0.5087; 
        %side stream inside radius [m] ~ 8"
        side_stream_radius =  0.2445; 
        %reservoir length, injection to production well [m]
        res_length;
        %Number of side streams
        n_streams = 1;
        %spacing between lateral streams [m]
        stream_spacing = 100;
        %angle of the lateral wells - horizontal by default [degrees]
        angle = 180; 
        %time [years]
        time_years = 30
        %surface average ambient air temp [C]
        T_surface_air_C = 10;
        %surface rock temp [C]
        T_surface_rock_C = 10;
        % earth's geothermal gradient [C/m], always positive
        dT_dz = 0.035;
        fluid = 'CO2';
        %fluid = 'Water';

        % CO2 properties
        % Specific heat ratio
        k = 1.289;

        pcrit = py.CoolProp.CoolProp.PropsSI('pcrit','CO2') %Pa
        Tcrit = py.CoolProp.CoolProp.PropsSI('Tcrit','CO2') - 273.15 %C

        %configurations [config == 1:4]
        config = 1;

        %cooling water system
        T_cooling_water = 10; %C
        P_cooling_water = 5e5; %Pa
        eta_cooling_water_pump = 0.75;

        %heat exchanger paramters
        U_dirty = 600; %W/m^2/K

        %splitter parameters

        S_ratio = 1 % Split ratio --> S = m1/m2; mtotal = m1+m2
        
        % toggle whether ideal S_ratio is calculated for config 3|| default
        % is off.
        find_opt_S_ratio_toggle = 0

        optimizationMode = 'MinimizeLCOE_Brownfield';
        %optimizationMode = 'MinimizeLCOE_Greenfield';
        %optimizationMode = 'MaximizePower';
        %silicaPrecip = 'PreventSilica';
        silicaPrecip = 'IgnoreSilica';
        coolingMode = 'Wet';
%         coolingMode = 'Dry';
        orcFluid = 'R245fa';
        %orcFluid = 'R600a';
        %thickness for porous heat depletion
        thickness = 100;
        %wellFieldType 
        %(Note: 5spot and 5spot_ManyN with N_5spot=1 are the same thing)
        %'Doublet','5spot','5spot_SharedNeighbor','5spot_ManyN'
        wellFieldType = 'Doublet';
        %N_5spot is only used with the '5spot_ManyN' wellfield type.
        %N_5spot is the sqrt of the number of 5spots laid adjacently
        N_5spot = 1;
        
        %O&M Fraction
        F_OM = 0.045;
        %Discount Rate
        discountRate = 0.096;
        %Financial Lifetime
        Lifetime = 25;
        %Capacity Factor
        CapacityFactor = 0.9;
        %Cost Year
        costYear = 2022;
        
        % 2ellCostType. Ideal values are the technological limit, Baseline
        % values are the current state of the art.
        % allowed: 'Ideal','Baseline'
        wellCostType = 'Baseline';
        % Success rate of wells
        SuccessRate_well = 0.95;
        
        %monitoring well diameter
        monitoringWellDiameter = 0.216;
        
        %porous reservoir simulation
        modelResPressureTransient = false;
        %porous reservoir simulation
        modelResTemperatureDepletion = true;
        
        %CPG params
        eta_cpg_turbine = 0.78;
        eta_cpg_pump = 0.9;
        eta_cpg_compressor = 0.78;

        %Turbine parameters config 2
        T_turb_out_design = 65; %C
        T_turb_out_min = 17; %C

        % Is phase change allowed in the turbine from supercritical
        % to liquid?
        PhaseChangeAllowed = false;
        
        % Overall heat transfer coefficient for orc HX
        HX_overallHeatTransferCoefficient = 500; %W/m^2-K

        % Pressure drop across a heat exchanger
        dP_hex = 0.5e5; % Pa (0.5 bar)

        % Refrigerant Parameters
        % Refrigerant
        R717 = 'R717'
        
        %ORC params
        dT_approach = 7;
        dT_orc_pinch = 5;
        eta_orc_pump = 0.9;
        eta_orc_turbine = 0.8;
        orcCycleType = 'Subcritical';
        %orcCycleType = 'Supercritical';
        orcModel = 'Simulation';
        %orcModel = 'Geophires';

        % hasSurfaceGatheringSystem (true/false). Sets if there is a system of
        % pipes connecting the surface injection and production wells.

        % The injection and production wells are nearby
        hasSurfaceGatheringSystem = false;
        
        %water params
        PumpDepth = 500;
        maxPump_dP = 10e6; %10 MPa max
        eta_pump = 0.75;
        
        %well parameters
        useWellboreHeatLoss = true;
        % These three should be inputs into this function
        k_rock = 2.96; %W/m-C
        C_rock = 1000; %J/kg-C
        rho_rock = 2650; %kg/m^3
        %Friction Factor
        epsilon = 55 * 1e-6; % 55 um

        %% Storage parameters
        %Storage Pressure
        P_store = 68e5
        %Storage Temperature
        T_store = 15;

        %% Process Pipe layout
        % Pipe length estimates
        % Main Process line (1 = turbine, 2 = domestic heating)
        surf_pipe_length_MP1 = 45 %m
        surf_pipe_length_MP2 = 35 %m
        % Main line start_up and storage
        surf_pipe_length_SUS = 140 %m
        % Cooling water 
        surf_pipe_length_CW = 105 %m

        % simulation in steady state 1 = yes, 0 = no (start-up conditions)
        steadystate = 1;

        %simulatorType 'genGEO' or 'GEOPHIRES'
        simulatorType = 'genGEO';

        %% FINANCIAL PARAMETERS
        PPI = PPIs.load_PPI();
   
    end
end