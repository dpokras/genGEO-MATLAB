classdef Surface_plant
    methods(Static)
        function P_surface_out = plant_dP(P_surface_in, T_surface_in, params)
            h_surface_in = CoolProp.PropsSI('HMASS', 'P', P_surface_in, 'T', T_surface_in+273.15, params.fluid);
            
            rho_surface = CoolProp.PropsSI('DMASS', 'P', P_surface_in, 'HMASS', h_surface_in, params.fluid);
            ff = FrictionFactor(params.well_radius, P_surface_in, h_surface_in, params.m_dot, params);
            if (params.hasSurfaceGatheringSystem == true)
                dP_surfacePipes_MP1 = ff * params.surf_pipe_length_MP1 / (params.well_radius*2)^5 * 8 * params.m_dot^2 / rho_surface / pi^2;
                dP_surfacePipes_MP2 = ff * params.surf_pipe_length_MP1 / (params.well_radius*2)^5 * 8 * params.m_dot^2 / rho_surface / pi^2;
                dP_surfacePipes = dP_surfacePipes_MP1 + dP_surfacePipes_MP2;
            else
                dP_surfacePipes = 0;
            end
             P_surface_out = P_surface_in - dP_surfacePipes;
        end

        function result = TurbWorkfunc(P_turb_in, T_turb_in, params)
            
            if params.config == 4
                W_turbine = 0;
                T_turb_out_irr = T_turb_in;
                P_turb_out_irr = P_turb_in;
                s_turb_in = CoolProp.PropsSI('SMASS', 'P', P_turb_in, 'T', T_turb_in+273.15, params.fluid);
            else
    
                %%% States
                %%State 1
                h_turb_in = CoolProp.PropsSI('HMASS', 'P', P_turb_in, 'T', T_turb_in+273.15, params.fluid);
                s_turb_in = CoolProp.PropsSI('SMASS', 'P', P_turb_in, 'T', T_turb_in+273.15, params.fluid);
    
                %%State 2
                if params.config == 2
                     P_turb_out_iter = 115e5;
                else
                    P_turb_out_iter = 74e5;
                end
                toggle = 0;
                go = 1;
                iter = 1;
    
                while go == 1
                
    %                 disp(strcat(['interation: ' num2str(iter)]));
                    
                    %isentropic output
                    h_turb_out_isen = CoolProp.PropsSI('HMASS', 'P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid);
                    
                    %efficiency losses
                    h_turb_out_irr = h_turb_in - params.eta_cpg_turbine * (h_turb_in - h_turb_out_isen);
                    s_turb_out_irr = CoolProp.PropsSI('SMASS', 'P', P_turb_out_iter, 'HMASS', h_turb_out_irr, params.fluid);
                    T_turb_out_irr = CoolProp.PropsSI('T', 'P', P_turb_out_iter, 'SMASS', s_turb_out_irr, params.fluid);
                    phase = CoolProp.PhaseSI('P', P_turb_out_iter, 'SMASS', s_turb_in, params.fluid);
                    
    %                 disp(strcat(['Turbine Temperature Out: ' num2str(T_turb_out_irr-273.15) 'C']));
    %                 disp(strcat(['Turbine Pressure Out: ' num2str(P_turb_out_iter/1e5) ' bar']));
    %                 disp(phase);
    
                    P_turb_out_iter = P_turb_out_iter - 0.1e5;
                    iter = iter + 1;
    
                    % Turbine outlet pressure is determined by the while loop above,
                    % and a 0.5 bar buffer is added to avoid phase change and 
                    % damage to the machinery.
                    %
                    % PhaseChangeAllowed toggle means if the turbine can handle
                    % phase change from supercritical fluid to liquid.
    
                    if toggle == 1
                        go = 0;
    
                    elseif params.config == 2
                        if T_turb_out_irr < params.T_turb_out_design + 273.15
                            P_turb_out_iter = P_turb_out_iter + 1.1e5;
                            toggle = 1;
                        end
                    else
                        if strcmp(string(phase), 'liquid') == 1 || contains(string(phase), 'twophase') == 1
        
                            P_turb_out_iter = P_turb_out_iter + 1.1e5;
                            toggle = 1;
                        end
                    end
                end
                P_turb_out_irr = P_turb_out_iter;
                if params.config == 3
                    m_dot_turb = params.m_dot * params.S_ratio;
                elseif params.config == 4
                    m_dot_turb = 0;
                else
                    m_dot_turb = params.m_dot;
                end

                W_turbine =  m_dot_turb * (h_turb_in - h_turb_out_irr);

                if (W_turbine < 0)
                    throw(MException('total_analytic_system_CO2:TurbinePowerNegative','Turbine Power is Negative'));
                end
            end
            result.W_turbine = W_turbine;
            result.T = T_turb_out_irr-273.15;
            result.P = P_turb_out_irr;
            result.s_turb_in = s_turb_in;
        end

        function result = Compfunc(params)

            T_comp_in = 30; %C
            T_comp_out = 50; %C
            P_comp_out = 70e5; %Pa

            h_comp_in = CoolProp.PropsSI('HMASS', 'P', params.P_store, 'T', T_comp_in + 273.15, params.fluid);
            h_comp_out_isen = CoolProp.PropsSI('HMASS', 'P', P_comp_out, 'T', T_comp_out + 273.15, params.fluid);
            h_comp_out_irr = h_comp_in - params.eta_cpg_compressor * (h_comp_in - h_comp_out_isen);
            
            % compressor should be designed to handle a flow at least 10%
            % of the total flow rate of the process.
            m_dot = params.m_dot * 0.1;
            W_compressor = m_dot * (h_comp_out_irr - h_comp_in);

            % calcualte the cooling HEX afterthe compressor
            P_hex_hot_in = P_comp_out; %Pa
            T_hex_hot_in = T_comp_out; %C
            T_hex_hot_out = 17; %C
            T_hex_cold_in = params.T_cooling_water; %C
            T_hex_cold_out = 40; %C
            
            h_hex_hot_in = CoolProp.PropsSI('HMASS', 'P', P_hex_hot_in, 'T', T_hex_hot_in + 273.15, params.fluid);
            h_hex_hot_out = CoolProp.PropsSI('HMASS', 'P', P_hex_hot_in - params.dP_hex, 'T', T_hex_hot_out + 273.15, params.fluid); %J/kg
            h_hex_cold_in = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex_cold_in + 273.15, 'water'); %J/kg
            h_hex_cold_out = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex_cold_out + 273.15, 'water'); %J/kg
            
            m_dot_water = (m_dot * (h_hex_hot_in-h_hex_hot_out) )/(h_hex_cold_out-h_hex_cold_in);

            % number of shells in series
            N = Surface_plant.HEX_shellcount(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out);
            LMTD = Surface_plant.LMTD(T_hex_hot_in, T_hex_hot_out, T_hex_cold_in, T_hex_cold_out, N, 0); %K
            Q = m_dot * (h_hex_hot_in - h_hex_hot_out); %W            
            A = Q / params.U_dirty / LMTD; %m^2

            result.A = A;
            result.Q = Q;
            result.P_tubeside = P_hex_hot_in;
            result.W_compressor = W_compressor;
            result.m_dot_water = m_dot_water;
        end
        

        %% Condenser Section
        % If phase is gas, then cooling and condensing gas to liquid
        % phase. If phase is liquid, then no condensation.

        function result = HEXfunc(P_turb_out_irr, T_turb_out_irr, T_well_out, P_well_out, params)
            
            if params.config == 1
                P_hex1_hot_in = P_turb_out_irr; %Pa
                T_hex1_hot_in = T_turb_out_irr; %C
                T_hex1_hot_out = CoolProp.PropsSI('T', 'P', P_hex1_hot_in - params.dP_hex, 'Q', 0, params.fluid) - 273.15; %J/kg
                T_hex1_cold_in = params.T_cooling_water; %C
                T_hex1_cold_out = T_hex1_hot_in - params.dT_approach; %C
                
                h_hex1_hot_in = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in, 'T', T_hex1_hot_in + 273.15, params.fluid);
                h_hex1_hot_satvap = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in - params.dP_hex, 'Q', 1, params.fluid); %J/kg
                h_hex1_hot_out = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in - params.dP_hex, 'Q', 0, params.fluid); %J/kg
                h_hex1_cold_in = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex1_cold_in + 273.15, 'water'); %J/kg
                h_hex1_cold_out = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex1_cold_out + 273.15, 'water'); %J/kg
                h_hex2_hot_out = 0;

                m_dot_water_CoolLoop = (params.m_dot * (h_hex1_hot_in-h_hex1_hot_out) )/(h_hex1_cold_out-h_hex1_cold_in);

                % proportion of desuperheating from total cooling process.
                cool_frac = (h_hex1_hot_in - h_hex1_hot_satvap) / (h_hex1_hot_in - h_hex1_hot_out);
                
                % number of shells in series
                N_1 = Surface_plant.HEX_shellcount(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out);
                LMTD_1 = Surface_plant.LMTD(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out, N_1, cool_frac); %K
                Q_hex_1 = params.m_dot * (h_hex1_hot_in - h_hex1_hot_out); %W
                Q_hex_2 = 0;
                Q_net = 0;

                result_cooling_tower = Surface_plant.CoolingTowerfunc(params.m_dot, P_hex1_hot_in, T_hex1_hot_in, h_hex1_hot_out, params);
                
                A_1 = Q_hex_1 / params.U_dirty / LMTD_1; %m^2
                A_2 = 0;
                N_2 = 0;

                T_CO2_out = T_hex1_hot_out;
                T_water_out = T_hex1_cold_out;
                P_CO2_out = P_hex1_hot_in - params.dP_hex;
                P_water_out = params.P_cooling_water - params.dP_hex;
                m_dot_water = 0;
                P_tubeside_1 = P_hex1_hot_in - params.dP_hex/2;
                P_tubeside_2 = 0;


            elseif params.config == 2 || params.config == 4 || params.S_ratio == 0

                if params.config == 2 
                    T_hex1_hot_in = T_turb_out_irr; %C
                    P_hex1_hot_in = P_turb_out_irr; %Pa
                elseif params.config == 4 || params.S_ratio == 0
                    T_hex1_hot_in = T_well_out; %C
                    P_hex1_hot_in = P_well_out; %Pa
                end
                
                T_hex1_hot_out = params.T_cooling_water + params.dT_approach; %C
                T_hex1_cold_in = params.T_cooling_water; %C
                T_hex1_cold_out = 50; %C

                % number of shells in series
                N_1 = Surface_plant.HEX_shellcount(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out);

                LMTD_1 = Surface_plant.LMTD(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out, N_1, 0); %K
                
                h_hex1_hot_in = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in, 'T', T_hex1_hot_in+273.15, params.fluid);
                h_hex1_hot_out = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in-params.dP_hex, 'T', T_hex1_hot_out+273.15-0.001, params.fluid); %J/kg
                h_hex1_cold_in = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex1_cold_in+273.15, 'water'); %J/kg
                h_hex1_cold_out = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex1_cold_out+273.15, 'water'); %J/kg
                h_hex2_hot_out = 0;
                Q_hex_1 = params.m_dot * (h_hex1_hot_in - h_hex1_hot_out); %W
                Q_hex_2 = 0;
                Q_net = Q_hex_1;

                result_cooling_tower.Q_desuperheating = 0; %W
                result_cooling_tower.Q_condensing = 0; %W
                result_cooling_tower.Q_cooling_tower = 0; %W
                result_cooling_tower.W_cooling_tower = 0; %W
                result_cooling_tower.dT_range = 0; %K
                result_cooling_tower.m_dot_water_makeup = 0; %kg/s
                
                A_1 = Q_hex_1 / params.U_dirty / LMTD_1; %m^2
                A_2 = 0;
                N_2 = 0;

                m_dot_water  = Q_hex_1 / (h_hex1_cold_out - h_hex1_cold_in);
                m_dot_water_CoolLoop = 0;
                if params.config == 4 || params.S_ratio == 0
                    P_CO2_out = 6.5e6; %Pa    
                    T_CO2_out = CoolProp.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex1_hot_out, params.fluid) -273.15; %J/kg
                else
                    P_CO2_out = P_hex1_hot_in - params.dP_hex;
                    T_CO2_out = T_hex1_hot_out;
                end
                T_water_out = T_hex1_cold_out;
                P_water_out = params.P_cooling_water - params.dP_hex;
                P_tubeside_1 = P_hex1_hot_in - params.dP_hex/2;
                P_tubeside_2 = 0;

            elseif params.config == 3

                m_dot_sco2_1 = params.m_dot*params.S_ratio;
                m_dot_sco2_2 = params.m_dot*(1-params.S_ratio);
                
                %HEX stream 1
                P_hex1_hot_in = P_turb_out_irr; %Pa
                T_hex1_hot_in = T_turb_out_irr; %C
                T_hex1_cold_in = params.T_cooling_water; %C
                T_hex1_cold_out = T_hex1_hot_in - params.dT_approach;
                
                h_hex1_hot_in = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in, 'T', T_hex1_hot_in + 273.15, params.fluid); %J/kg
                h_hex1_cold_in = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water, 'T', T_hex1_cold_in + 273.15, 'water'); %J/kg
                h_hex1_cold_out = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water - params.dP_hex, 'T', T_hex1_cold_out + 273.15, 'water'); %J/kg

                %HEX stream 2
                P_hex2_hot_in = P_well_out; %Pa
                T_hex2_cold_in = T_hex1_cold_out;
                T_hex2_hot_in = T_well_out; %C
                T_hex2_hot_out = T_hex2_cold_in + params.dT_approach; %C
                T_hex2_cold_out = 50; %C

                h_hex2_hot_in = CoolProp.PropsSI('HMASS', 'P', P_hex2_hot_in, 'T', T_hex2_hot_in + 273.15, params.fluid); %J/kg
                h_hex2_hot_out = CoolProp.PropsSI('HMASS', 'P', P_hex2_hot_in - params.dP_hex, 'T', T_hex2_hot_out + 273.15, params.fluid); %J/kg
                h_hex2_cold_out = CoolProp.PropsSI('HMASS', 'P', params.P_cooling_water - 2*params.dP_hex, 'T', T_hex2_cold_out + 273.15, 'water'); %J/kg
                
                % Combined streams
                h_hex1_hot_satvap = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in, 'Q', 1, params.fluid); %J/kg
                h_hex_hot_satliq = CoolProp.PropsSI('HMASS', 'P', P_hex1_hot_in-params.dP_hex, 'Q', 0, params.fluid); %J/kg

                P_CO2_out = P_hex1_hot_in - params.dP_hex;
                T_hot_satliq = CoolProp.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex_hot_satliq, params.fluid) - 273.15; %C

                if params.find_opt_S_ratio_toggle == 0

                    h_hex1_hot_out = m_dot_sco2_2/m_dot_sco2_1 * (h_hex_hot_satliq - h_hex2_hot_out) + h_hex_hot_satliq;
                    if m_dot_sco2_1 > 0
                        T_hex1_hot_out = CoolProp.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex1_hot_out, params.fluid) - 273.15; %C
                    end
                    T_CO2_out = T_hot_satliq;
                    
                    %%%% if below statement is true then this is incorrect and will be recalculated. h_hex1_hot_out is calculated as the required enthaply to
                    %%%% reach saturated liquid state when mixed with stream 2 output. In the case of a low S ratio --> than more subcooled liquid is in the
                    %%%% mix so a hoter stream 1 is technically required. Therefore no more cooling is need and is fully accompliushed with stream 2.
                    %%%% Therefore the final enthalpy is simply calculated after mixing the would be output streams.
                    if h_hex1_hot_out > h_hex1_hot_in
                        cooling_required = 0;
                    else
                        cooling_required = 1;
                    end

                elseif params.find_opt_S_ratio_toggle == 1

                    cooling_required = 0;
                    
                    params.S_ratio = (h_hex_hot_satliq - h_hex2_hot_out)/(h_hex1_hot_in - h_hex2_hot_out);

                    m_dot_sco2_1 = params.m_dot*params.S_ratio;
                    m_dot_sco2_2 = params.m_dot*(1-params.S_ratio);
                    
                    T_hex2_cold_in = T_hex1_cold_in;
                    h_hex2_cold_in = h_hex1_cold_in;
                    h_hex1_hot_out = h_hex1_hot_in;
                    m_dot_water = m_dot_sco2_2 * (h_hex2_hot_in-h_hex2_hot_out)/(h_hex2_cold_out-h_hex2_cold_in);
                    m_dot_water_CoolLoop = 0;
                end
                
                if cooling_required == 0
                    h_hex_hot_combined_out = params.S_ratio * h_hex1_hot_in + (1 - params.S_ratio)*h_hex2_hot_out;
                    T_CO2_out = CoolProp.PropsSI('T', 'P', P_CO2_out, 'HMASS', h_hex_hot_combined_out, params.fluid) - 273.15; %C
                end

                if params.find_opt_S_ratio_toggle == 1
                    
                elseif cooling_required == 0

                    h_hex2_cold_in = h_hex1_cold_in;
                    T_hex2_cold_in = T_hex1_cold_in;
                    %m_dot of water being sold for waste heat
                    m_dot_water = m_dot_sco2_2 * (h_hex2_hot_in-h_hex2_hot_out)/(h_hex2_cold_out-h_hex2_cold_in);
                    m_dot_water_CoolLoop = 0;

                elseif cooling_required == 1

                    %m_dot_water is water being sold for waste heat                        
                    m_dot_water = (m_dot_sco2_2 * (h_hex2_hot_in-h_hex2_hot_out))/(h_hex2_cold_out-h_hex1_cold_in);
                    m_dot_water_CoolLoop = m_dot_sco2_1 * (h_hex1_hot_in-h_hex1_hot_out)/(h_hex1_cold_out-h_hex1_cold_in);
                    h_hex2_cold_in = h_hex1_cold_out;
                    T_hex2_cold_in = CoolProp.PropsSI('T', 'P', params.P_cooling_water - params.dP_hex, 'HMASS', h_hex2_cold_in, 'water') - 273.15; %C
                    
                    cool_frac_1 = (h_hex1_hot_in - h_hex1_hot_satvap) / (h_hex1_hot_in - h_hex1_hot_out);

                end
                    
                if params.find_opt_S_ratio_toggle == 0 && cooling_required == 1
                    %%%% Heat Excahnger 1 (small) 
                    % number of shells in series
                    N_1 = Surface_plant.HEX_shellcount(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out);

                    LMTD_1 = Surface_plant.LMTD(T_hex1_hot_in, T_hex1_hot_out, T_hex1_cold_in, T_hex1_cold_out, N_1, cool_frac_1); %K
                    Q_hex_1 = m_dot_sco2_1 * (h_hex1_hot_in - h_hex1_hot_out); %W 
                    A_1 = Q_hex_1 / params.U_dirty / LMTD_1; %m^2
                    P_water_out = params.P_cooling_water - params.dP_hex;
                    P_tubeside_1 = P_hex1_hot_in - params.dP_hex/2;

                    result_cooling_tower = Surface_plant.CoolingTowerfunc(m_dot_sco2_1, P_hex1_hot_in, T_hex1_hot_in, h_hex1_hot_out, params);
                
                else
                    N_1 = 0;
                    Q_hex_1 = 0; %W
                    result_cooling_tower.Q_desuperheating = 0; %W
                    result_cooling_tower.Q_condensing = 0; %W
                    result_cooling_tower.Q_cooling_tower = 0; %W
                    result_cooling_tower.W_cooling_tower = 0; %W
                    result_cooling_tower.dT_range = 0; %K
                    result_cooling_tower.m_dot_water_makeup = 0; %kg/s

                    A_1 = 0; %m^2
                    P_water_out = params.P_cooling_water;
                    P_tubeside_1 = 0;
                end
                
                %%%% Heat Excahnger 2 (large)
                % number of shells in series
                N_2 = Surface_plant.HEX_shellcount(T_hex2_hot_in, T_hex2_hot_out, T_hex2_cold_in, T_hex2_cold_out);

                LMTD_2 = Surface_plant.LMTD(T_hex2_hot_in, T_hex2_hot_out, T_hex2_cold_in, T_hex2_cold_out, N_2, 0); %K
                Q_hex_2 = m_dot_sco2_2 * (h_hex2_hot_in - h_hex2_hot_out); %W
                A_2 = Q_hex_2 / params.U_dirty / LMTD_2; %m^2

                Q_net = m_dot_water * ((h_hex1_cold_out - h_hex1_cold_in) + (h_hex2_cold_out - h_hex2_cold_in));

                T_water_out = T_hex2_cold_out;
                P_CO2_out = P_hex1_hot_in - params.dP_hex;
                P_water_out = P_water_out - params.dP_hex;
                P_tubeside_2 = P_hex2_hot_in - params.dP_hex/2;
            end    

            result.Q_hex_1 = Q_hex_1;
            result.Q_hex_2 = Q_hex_2;
            result.result_cooling_tower = result_cooling_tower;
      
            %%% increasing the  the surface area by 50% to account for
            %%% fluctuations in cooling water inlet temperature and
            %%% thus reduced heat exchange efficiency
            
            result.Q_net = Q_net;
            result.A_1 = A_1 * 1.5;
            result.A_2 = A_2 * 1.5;
            result.N_1 = N_1;
            result.N_2 = N_2;
            result.m_dot_water = m_dot_water;
            result.m_dot_water_CoolLoop = m_dot_water_CoolLoop + result_cooling_tower.m_dot_water_makeup;
            result.T_water_out = T_water_out;
            result.T_CO2_out = T_CO2_out;
            result.P_water_out = P_water_out;
            result.P_CO2_out = P_CO2_out;
            result.P_tubeside_1 = P_tubeside_1;
            result.P_tubeside_2 = P_tubeside_2;
            result.h_hex1_hot_out = h_hex1_hot_out;
            result.h_hex2_hot_out = h_hex2_hot_out;
            result.S_ratio = params.S_ratio;

        end

        function LMTD = LMTD(T_hot_in, T_hot_out, T_cold_in, T_cold_out, N, cool_frac)
            % N is the number of shells in series
            dT_in = abs(T_hot_in-T_cold_out); %K
            dT_out = abs(T_hot_out-T_cold_in); %K
            LMTD = (dT_in-dT_out)/log(dT_in/dT_out);

            R = (T_hot_in - T_hot_out)/(T_cold_out-T_cold_in);
            P = (T_cold_out - T_cold_in)/(T_hot_in - T_cold_in);
            
            S = (R^2 + 1)^0.5 / (R - 1);
            W = ((1 - P*R)/(1 - P))^(1/N);
            F_T = (S*log(W))/ (log(( 1 + W - S + S*W) /( 1 + W + S - S*W)));

            if T_hot_in == T_hot_out
            elseif cool_frac > 0
                LMTD = LMTD*(cool_frac*F_T+(1-cool_frac)); 
            else
                LMTD = LMTD*F_T;
            end
        end
        function N = HEX_shellcount(T_hot_in, T_hot_out, T_cold_in, T_cold_out)
            % N is the number of shells in series. Based on graphical
            % technique
            m_cold = T_cold_out - T_cold_in;
            m_hot = T_hot_in - T_hot_out;
            
            y = T_hot_out;
            N = 0;
            x = 0;
            while x < 0.95
                x = (y - T_cold_in) / m_cold;
                y = m_hot * x + T_hot_out;
                N = N + 1;
            end
        end

        function result = RFHEXfunc(m_dot_cooling_water, T_cooling_water, P_cooling_water, params)
            
            % Absorber
            P_hex1_in_R717 = 1.2e5; %Pa
            T_hex1_in_R717 = -30; %C
            T_hex1_out_R717 = -10;

            P_hex1_in_water = P_cooling_water;
            T_hex1_in_water = T_cooling_water; %C
            T_hex1_out_water = params.T_cooling_water; %C

            h_hex1_in_R717 = CoolProp.PropsSI('HMASS', 'P', P_hex1_in_R717, 'T', T_hex1_in_R717 + 273.15, params.R717);
            h_hex1_out_R717 = CoolProp.PropsSI('HMASS', 'P', P_hex1_in_R717 - params.dP_hex, 'T', T_hex1_out_R717 + 273.15, params.R717); %J/kg
            h_hex1_in_water = CoolProp.PropsSI('HMASS', 'P', P_hex1_in_water, 'T', T_hex1_in_water + 273.15, 'water'); %J/kg
            h_hex1_out_water = CoolProp.PropsSI('HMASS', 'P', P_hex1_in_water - params.dP_hex, 'T', T_hex1_out_water + 273.15, 'water'); %J/kg
            
            m_dot_R717 = (m_dot_cooling_water * (h_hex1_in_water - h_hex1_out_water))/(h_hex1_out_R717 - h_hex1_in_R717);

            %Evaporator
            P_hex2_in_R717 = 20e5; %Pa
            T_hex2_in_R717 = 252; %C
            T_hex2_out_R717 = 50; %C

            P_hex2_air = 1.013e5; %Pa
            T_hex2_in_air = 25; %C
            T_hex2_out_air = 40; %C

            h_hex2_in_R717 = CoolProp.PropsSI('HMASS', 'P', P_hex2_in_R717, 'T', T_hex2_in_R717 + 273.15, params.R717);
            h_hex2_out_R717 = CoolProp.PropsSI('HMASS', 'P', P_hex2_in_R717 - params.dP_hex, 'T', T_hex2_out_R717 + 273.15, params.R717); %J/kg
            h_hex2_in_air = CoolProp.PropsSI('HMASS', 'P', P_hex2_air, 'T', T_hex2_in_air + 273.15, 'air'); %J/kg
            h_hex2_out_air = CoolProp.PropsSI('HMASS', 'P', P_hex2_air, 'T', T_hex2_out_air + 273.15, 'air'); %J/kg

            % proportion of desuperheating from total cooling process.
            m_dot_air = (m_dot_R717 * (h_hex2_in_R717 - h_hex2_out_R717))/(h_hex2_out_air - h_hex2_in_air);

            result.m_dot_air = m_dot_air;
            result.m_dot_R717 = m_dot_R717;

            %% We will not be including refrigeration because not enough
            % time to calculate the refrigeration cycle. Come back when
            % have more time. Out of scpe for thesis.
        end

        
        function result = CoolingTowerfunc(m_dot, P_cond_in, T_cond_in, h_cond_out, params)
            
            %heat rejection
            h_cond_in = CoolProp.PropsSI('HMASS', 'P', P_cond_in, 'T', T_cond_in + 273.15, params.fluid);

            if P_cond_in < params.pcrit
                h_satVapor = CoolProp.PropsSI('HMASS', 'P', P_cond_in, 'Q', 1, params.fluid);
            else
                h_satVapor = 1e7;
            end

            T_cond_out = CoolProp.PropsSI('T', 'P', P_cond_in - params.dP_hex, 'HMASS', h_cond_out, params.fluid) - 273.15;

            if h_cond_out > h_satVapor || P_cond_in > params.pcrit
                % only desuperheating needed, no condensation required
                Q_cooler_part = m_dot * (h_cond_in - h_cond_out);
                Q_condenser_part = 0;
                dT_range = T_cond_in - T_cond_out;
                m_dot_makeup = 0;

            elseif (h_cond_in > h_satVapor)
                %desuperheating needed
                Q_cooler_part = m_dot * (h_cond_in - h_satVapor);
                Q_condenser_part = m_dot * (h_satVapor - h_cond_out);
                dT_range = T_cond_in - T_cond_out;
                m_dot_makeup = (Q_cooler_part + Q_condenser_part)/(h_satVapor - h_cond_out);
            else
                %no desuperheating
                Q_cooler_part = 0;
                Q_condenser_part = m_dot * (h_cond_in - h_cond_out);
                dT_range = 0;
                m_dot_makeup = (Q_cooler_part+Q_condenser_part)/(h_cond_in - h_cond_out);
            end
            
            [f_cooling, f_condensing] = ParasiticPowerFraction_CoolingTower(params.T_surface_air_C, params.dT_approach, dT_range, params.coolingMode);
            W_cooler_part = f_cooling * Q_cooler_part;
            W_condenser_part = f_condensing * Q_condenser_part;
            Q_cooling_tower = Q_cooler_part + Q_condenser_part; 
            W_cooling_tower = W_cooler_part + W_condenser_part;
            
            result.Q_desuperheating = Q_cooler_part;
            result.Q_condensing = Q_condenser_part;
            result.Q_cooling_tower = Q_cooling_tower;
            result.W_cooling_tower = W_cooling_tower;
            result.dT_range = dT_range;
            result.m_dot_water_makeup = m_dot_makeup;

        end

        function result = PumpWorkfunc(dP, m_dot, P, T, fluid, params)

            % calculate net pump power value
            if strcmp(string(fluid), 'water') == 1
                rho = CoolProp.PropsSI('DMASS', 'P', P, 'T', T + 273.15, fluid); %kg/m^3
            elseif strcmp(string(fluid), 'CO2') == 1
                rho = CoolProp.PropsSI('DMASS', 'P', P, 'T', T + 273.15, fluid); %kg/m^3
            end
            
            Q = m_dot / rho * 15850; % gallons per minute;
            H = (dP/1e5) * 10.197 * 3.281; %ft

            % Electric motor calculations for pump
            mu_P = -0.316 + 0.24015*(log(Q)) - 0.01199*(log(Q))^2;
            P_B = Q*H*(rho/119.8)/33000/mu_P; %Hp
            mu_M = 0.80 + 0.0319*(log(P_B)) - 0.00182*(log(P_B))^2;
            P_c = P_B/mu_M; % Hp

            W_pump = P_c * 745.7;
            
            result.W_pump = W_pump;
            result.dP = dP;
            result.m_dot = m_dot;
            result.rho = rho;

        end

        function result = tank(injWell, reservoir, prodWell, params)

            m1 = pi * ((repmat(injWell.wellRadius.', params.n_streams,1)).').^2 .* injWell.Density(2:end,:) .* injWell.dL;
            m2 = pi * ((repmat(reservoir.wellRadius.', params.n_streams,1)).').^2 .* reservoir.Density(2:end,:) .* reservoir.dL.*ones(1, params.n_streams);
            m3 = pi * ((repmat(prodWell.wellRadius.', params.n_streams,1)).').^2 .* prodWell.Density(2:end,:) .* prodWell.dL;

            m_total = sum(m1(:,end)) + sum(m2,'all') + sum(m3(:,end)); % kg

            % Rho at storage condition: P above critical to avoid phase change and rapid decompression         
            rho_store = CoolProp.PropsSI('DMASS','P', params.P_store, 'T', params.T_store+273.15, params.fluid); %kg/m^3
            
            %Store 10% for contingency
            V_total_surface = m_total / rho_store; %m^3
            V_store = V_total_surface * 0.1;

            result.m_total = m_total;
            result.V_total_surface = V_total_surface;
            result.V_store = V_store;
        end

        function W_net = NetWorkfunc(W_turbine, W_cooling_tower, W_storage_comp, result_pumps, params)
            
            W_cooling_pump = result_pumps.result_cooling_pump.W_pump;
            W_SU_pump = result_pumps.result_SU_pump.W_pump;
            W_filling_pump = result_pumps.result_filling_pump.W_pump;

            if params.steadystate == 1
            % calculate net power values
                W_net = W_turbine - W_cooling_pump - W_cooling_tower;
            else
                W_net = W_turbine - W_cooling_pump - W_cooling_tower - W_storage_comp - W_SU_pump - W_filling_pump;
            end

        end
    end
end