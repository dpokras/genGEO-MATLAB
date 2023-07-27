classdef CapitalCost_SurfacePlantfunc
    methods(Static)
        function C_BM = HEX(area, P, params)
            %% STEADY STATE COSTS

            %% USD to EURO Conversion (May 2023)
            % 1 USD = 0.92 EURO

            %% Heat exchangers
            % Description: fixed head type Shell & tube HEX. TEMA type DEM Carbon steel except 
            % tube-side Cr-Mo steel for corrosion resistance. sCO2 is tube-side; 
            % H2O is shell-side.

            % REFERENCE: Product and Process Design Principles: Synthesis,
            % Analysis, and Evaluation 4th Ed.
            
            % Heat exchanger 1
            
            % reference CEPCI; 2013 (average)
            PPI_ref = PPIs.return_PPI('PPI_HX', 2013, params); %2013 330.167
            PPI_current = PPIs.return_PPI('PPI_HX', params.costYear, params); %2022 494.764
                
            % inputs
            area = area .* 10.764; % ft^2

            % free-of-board costs
            C_B = exp(11.4185 - 0.9228.*log(area) + 0.09861.*(log(area)).^2); %USD
            C_B(isinf(C_B)) = 0;
            
            % material factor -- S/T-side = carbon steel/Cr-Mo steel 
            a = 1.55;
            b = 0.05;
            F_M = a + (area./100).^b;
            
            % length factor L = 6069mm
            F_L = 1.0;
            
            % tubeside pressure factor - cost correlation (P > 10.342 bar)
            % REFERENCE: CAMARAZA-MEDINA et al. (2021) https://doi.org/10.52292/j.laar.2021.713
            if P > 10.342e5 % Pa
                %Inner diameter is assumed to be 650mm.
                ID = 0.65; %m
                F_P = 1 + ((P./1e5)./10.342 - 1).*(0.0035 - 0.022 .* (ID-0.3048)); % P in Pa, ID in m
            else
                F_P = 1.0;
            end
            
            C_P = F_P.*F_M.*F_L.*C_B; %USD
            C_P = C_P .* (PPI_current./PPI_ref); %USD

            % Guthrie Bare-Module Factor for STHEs
            F_BM = 3.17;
            C_BM = C_P * F_BM;


        end
        
        function C_BM = pump(input, params)
            %% Pump Cost Estimation
        
            % reference CEPCI; 2013 (average)
            PPI_ref = PPIs.return_PPI('PPI_Pump&Comp', 2013, params); %2013 140.167
            PPI_current = PPIs.return_PPI('PPI_Pump&Comp', params.costYear, params); %2022 192.573
            
            Q = input.m_dot / input.rho * 15850; % gallons per minute;
            H = (input.dP/1e5) * 10.197 * 3.281; %ft
            
            S = Q*(H)^0.5; % Size factor
            
            C_B_pump = exp(12.1656 - 1.1448*log(S)+0.0862*(log(S))^2); %USD
            
            % pump type factor; 1 stage, 1800rpm, 250 - 5000gpm range, 50 - 500ft head range, 250Hp max motor power 
            F_T_pump = 2.0;
            
            % material factor; cast iron
            F_M_pump = 1.0;
            
            % Electric motor calculations for pump
            mu_P = -0.316 + 0.24015*(log(Q)) - 0.01199*(log(Q))^2;
            P_B = Q*H*(input.rho/119.8)/33000/mu_P; %Hp
            mu_M = 0.80 + 0.0319*(log(P_B)) - 0.00182*(log(P_B))^2;
            P_c = P_B/mu_M; % Hp
            
            C_B_driver = exp(5.9332 + 0.16829*log(P_c) - 0.110056*(log(P_c))^2 + ...
                0.071413 * (log(P_c))^3 - 0.0063788 * (log(P_c))^4); %USD
            
            F_T_driver = 1.3; %1800rpm enclosed, fan-cooled, 1 to  250 Hp

            C_P = C_B_pump*F_T_pump*F_M_pump + C_B_driver*F_T_driver; %USD
            C_P = C_P * (PPI_current / PPI_ref); %USD

            % Guthrie Bare-Module Factor for pumps
            F_BM = 3.30;
            C_BM = C_P * F_BM;

        end

        function C_BM = compressor(input, params)
            %% Compressor Cost Estimation
            
            % reference CEPCI; 2013 (average)
            PPI_ref = PPIs.return_PPI('PPI_Pump&Comp', 2013, params); %2013 140.167
            PPI_current = PPIs.return_PPI('PPI_Pump&Comp', params.costYear, params); %2022 192.573

            P_c = input.W_compressor / 745.7; % Hp
            C_B = exp(9.1553 + 0.63*log(P_c)); %USD
            
            F_D = 1.0;
            F_M = 2.5;
            
            C_P = C_B * F_M * F_D; %USD
            C_P = C_P * (PPI_current / PPI_ref); %USD

            % Guthrie Bare-Module Factor for gas compressors
            F_BM = 2.15;
            C_BM = C_P * F_BM;

        end
        function C_BM = turbine(input, params)

            % reference CEPCI; 2002 (average)
            PPI_ref = PPIs.return_PPI('PPI_T-G', 2003, params); %2003 154.033
            PPI_current = PPIs.return_PPI('PPI_T-G', params.costYear, params); %2022 245.235            
            % Calculate:
            % C_T_G
            %Regular fluid
            S_T_fluid = 1.20; %CO2
            C_P = 0.67 * (S_T_fluid*2830*(input.W_turbine/1e3)^0.745 + 3680*(input.W_turbine/1e3)^0.617);
            
            C_P = C_P * (PPI_current / PPI_ref); %USD

            % Guthrie Bare-Module Factor for gas_driven turbines
            F_BM = 1.5;
            C_BM = C_P * F_BM;

        end
        function C_BM = cooling_tower(input, params)

            % reference CEPCI; 2019 (average)
            PPI_ref = PPIs.return_PPI('PPI_Cool', 2019, params); %2013 131.650
            PPI_current = PPIs.return_PPI('PPI_Cool', params.costYear, params); %2022 272.893

%             % Alternative calculation: Air-cooled fin-fan HEX
%             % REFERENCE: Product and Process Design Principles
%             
%             area = area * 10.764; % ft^2
%             C_P = 2835*area^0.40; %USD

            % Tower Design Coefficient - CO2: 1.20|H2O: 0.252
            TDC = 0.252;
            C_P = CapitalCost_CoolingTower(input.Q_desuperheating, ...
                input.Q_condensing, TDC, input.dT_range, params);

            C_P = C_P * (PPI_current / PPI_ref); %USD

            % Guthrie Bare-Module Factor for air-fin coolers
            F_BM = 2.17;
            C_BM = C_P * F_BM; %USD

        end

        function C_BM = tank(input, params)
            
            % REFERENCE: Product and Process Design Principles: Synthesis,
            % Analysis, and Evaluation 4th Ed.

            % reference CEPCI; 2013 (average)
            PPI_ref = PPIs.return_PPI('PPI_Tank', 2013, params); %2013 100.00
            PPI_current = PPIs.return_PPI('PPI_Tank', params.costYear, params); %2022 232.443 
            
            %V [m^3], P [Pa]
            % 20% of volume reserved for vapour --> V * 1.2
            V = input.V_store * 1.2; %m^3
            P = params.P_store/6895; %psi

            P_d = exp(0.60608 + 0.91615*log(P)+0.0015655*(log(P)^2));

            S = 15000; %psi
            E = 0.85;
%             t_s = 0.0063; %m
%             CA = 0.00315/39.3701; %inches

            % D:L ratio = 1:3
            % max diameter = 4m
            N = 8; %number of vessels
            D = ((V/N)*4/3/pi)^(1/3); %m
            D = D * 3.28084; %ft
            L = D * 3; %ft

            t_s = (P_d * (D*12) / (2*S*E - 1.2*P_d)); %inches
            
            rho_V = 0.284; %lb/in^3
            W = pi * ((D*12) + t_s)*((L*12) + 0.8*(D*12))*t_s*rho_V; %lb

            C_V = exp(5.6336 + 0.4599*log(W)+0.00582*(log(W))^2);
            C_V = C_V * N;
            
            % valid for D/L values between: 3<D<21; 12<L<40
%             C_PL = 410*((D)^0.73960)*((L)^0.70684);
            C_PL = 2275*D^0.2094;

            % Guthrie Bare-Module Factors for vertical pressure vessels
            F_BM = 3.05;

            %carbon steel
            F_M = 1.0;

            C_P = C_V * F_M + C_PL; %USD
            C_BM = C_P * F_BM; %USD
            C_BM = C_BM * (PPI_current / PPI_ref); %USD
            
        end

        function C_BM = cooling_water(input, params)
            
            CEPCI_current = PPIs.return_PPI('CEPCI', params.costYear, params); %2022 816

            rho = CoolProp.PropsSI('DMASS', 'P', P_pump_in, 'T', T_pump_in); %kg/m^3
            q = m_dot / rho; % m^3/s
            % 0.01 <  q < 10 [m^3/s]
            % REFERENCE: 

            a = 0.00007 + 2.5e-5 / q;
            b = 0.003;

            C_BM_unit = a * CEPCI_current + b; %USD/m^3
            
            % 2.88e7 s in 8000h (1 year)
            C_BM = C_BM_unit*q*2.88e7; %USD/year
        end

        function C_BM = refrigeration(input, params)

            CEPCI_current = PPIs.return_PPI('CEPCI', params.costYear, params); %2022 816

            % 0<T<300 [K], 0 <  Q_c < 1000 [kJ/s]
            % REFERENCE: 

            a = 0.5 * (Q_c/1e3)^(-0.9) / (T^3);
            b = 1.1e6 / T^5;
            
            C_BM_unit = a * CEPCI_current + b; %USD/kJ

            C_BM = C_BM_unit*(Q_c/1e3)*2.88e7; %USD/year
        end

        function C_TBM = external_costs(C_TBM)
            
            % site costs 10-20% total bare-module costs (TBM) for grass-roots plants
            % 4-6% for additon to integrated complex.
            % Assumed 10% of TBM.
            C_site = 0.1*C_TBM;
            
            % non-process buildings cost are: (1) 10% of TBM if process
            % equipment is housed (assumed to be the case); (2) 20% if 
            % un-housed and grassroots plant; (3) addtion to existing
            % plant.
            C_buildings = 0.2 * C_TBM;

            % offsite facilities are 5% of TBM.
            C_offsite = 0.05 * C_TBM;

            C_TBM = C_TBM + C_site + C_buildings + C_offsite;
            
        end
        function C_TCI = TCI(C_TBM)
            % 1.18 factor includes 15% contingency and 3% contractor fee
            C_TPI = 1.18*C_TBM;
            
            % Investment Site Factor based on location
            F_ISF = 1.20;
            C_TPI = C_TPI * F_ISF;
            
            % Working captial is 15% of Total Capital Investment or 17.6%
            % of total permanant investment.
            C_WC = 0.176 * C_TPI;
            
            % Total Captial Investment
            C_TCI = C_TPI + C_WC;
        end
    end
end
