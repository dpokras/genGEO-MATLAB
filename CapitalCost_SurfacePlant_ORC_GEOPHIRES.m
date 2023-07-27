function [C_plant] = CapitalCost_SurfacePlant_ORC_GEOPHIRES(T_in_C, ElectricityProduced_MW, cycleType, params)
    % Electricity in MW
    % T in C
    % Cost in dollars

    % GEOPHIRES Relations
    if (strcmp(cycleType,'SubcriticalORC'))
        MaxProducedTemperature = T_in_C;
        if (MaxProducedTemperature < 150)
            C3 = -1.458333E-3;
            C2 = 7.6875E-1;
            C1 = -1.347917E2; 
            C0 = 1.0075E4;
            CCAPP1 = C3*MaxProducedTemperature^3 + C2*MaxProducedTemperature^2 + C1*MaxProducedTemperature + C0;
        else
            CCAPP1 = 2231 - 2*(MaxProducedTemperature-150);
        end
        Cplantcorrelation = CCAPP1*(ElectricityProduced_MW/15)^(-0.06)*ElectricityProduced_MW*1000./1E6;
    elseif (strcmp(cycleType,'SupercriticalORC'))
        MaxProducedTemperature = T_in_C;
        if (MaxProducedTemperature < 150)
            C3 = -1.458333E-3;
            C2 = 7.6875E-1;
            C1 = -1.347917E2;
            C0 = 1.0075E4;
            CCAPP1 = C3*MaxProducedTemperature^3 + C2*MaxProducedTemperature^2 + C1*MaxProducedTemperature + C0;
        else
            CCAPP1 = 2231 - 2*(MaxProducedTemperature-150);
        end
        %factor 1.1 to make supercritical 10% more expansive than subcritical
        Cplantcorrelation = 1.1*CCAPP1*(ElectricityProduced_MW/15)^(-0.06)*ElectricityProduced_MW*1000./1E6;
    elseif (strcmp(cycleType,'SingleFlash'))
        if (ElectricityProduced_MW<10)
            C2 = 4.8472E-2; 
            C1 = -35.2186;
            C0 = 8.4474E3;
            D2 = 4.0604E-2; 
            D1 = -29.3817;
            D0 = 6.9911E3;
            PLL = 5;
            PRL = 10;
        elseif (ElectricityProduced_MW<25)
            C2 = 4.0604E-2; 
            C1 = -29.3817;
            C0 = 6.9911E3;	  
            D2 = 3.2773E-2; 
            D1 = -23.5519;
            D0 = 5.5263E3;        
            PLL = 10;
            PRL = 25;
        elseif (ElectricityProduced_MW<50)
            C2 = 3.2773E-2; 
            C1 = -23.5519;
            C0 = 5.5263E3;
            D2 = 3.4716E-2; 
            D1 = -23.8139;
            D0 = 5.1787E3;	          
            PLL = 25;
            PRL = 50;
        elseif (ElectricityProduced_MW<75)
            C2 = 3.4716E-2; 
            C1 = -23.8139;
            C0 = 5.1787E3;
            D2 = 3.5271E-2; 
            D1 = -24.3962;
            D0 = 5.1972E3;
            PLL = 50;
            PRL = 75;
        else
            C2 = 3.5271E-2; 
            C1 = -24.3962;
            C0 = 5.1972E3;
            D2 = 3.3908E-2;
            D1 = -23.4890;
            D0 = 5.0238E3;
            PLL = 75;
            PRL = 100;
        end
        maxProdTemp = T_in_C;
        CCAPPLL = C2*maxProdTemp^2 + C1*maxProdTemp + C0;
        CCAPPRL = D2*maxProdTemp^2 + D1*maxProdTemp + D0;
        b = log(CCAPPRL/CCAPPLL)/log(PRL/PLL);
        a = CCAPPRL/PRL^b;
        %factor 0.75 to make double flash 25% more expansive than single flash
        Cplantcorrelation = 0.8*a*(ElectricityProduced_MW^b)*ElectricityProduced_MW*1000./1E6;
        
    elseif (strcmp(cycleType,'DoubleFlash'))
        if (ElectricityProduced_MW<10)
            C2 = 4.8472E-2; 
            C1 = -35.2186;
            C0 = 8.4474E3;
            D2 = 4.0604E-2; 
            D1 = -29.3817;
            D0 = 6.9911E3;
            PLL = 5;
            PRL = 10;
        elseif (ElectricityProduced_MW<25)
            C2 = 4.0604E-2; 
            C1 = -29.3817;
            C0 = 6.9911E3;	  
            D2 = 3.2773E-2; 
            D1 = -23.5519;
            D0 = 5.5263E3;        
            PLL = 10;
            PRL = 25;
        elseif (ElectricityProduced_MW<50)
            C2 = 3.2773E-2; 
            C1 = -23.5519;
            C0 = 5.5263E3;
            D2 = 3.4716E-2; 
            D1 = -23.8139;
            D0 = 5.1787E3;	          
            PLL = 25;
            PRL = 50;
        elseif (ElectricityProduced_MW<75)
            C2 = 3.4716E-2; 
            C1 = -23.8139;
            C0 = 5.1787E3;	
            D2 = 3.5271E-2; 
            D1 = -24.3962;
            D0 = 5.1972E3;          
            PLL = 50;
            PRL = 75;
        else
            C2 = 3.5271E-2; 
            C1 = -24.3962;
            C0 = 5.1972E3;	
            D2 = 3.3908E-2; 
            D1 = -23.4890;
            D0 = 5.0238E3;          
            PLL = 75;
            PRL = 100;
        end
        maxProdTemp = T_in_C;
        CCAPPLL = C2*maxProdTemp^2 + C1*maxProdTemp + C0;
        CCAPPRL = D2*maxProdTemp^2 + D1*maxProdTemp + D0;
        b = log(CCAPPRL/CCAPPLL)/log(PRL/PLL);
        a = CCAPPRL/PRL^b;
        Cplantcorrelation = a*(ElectricityProduced_MW^b)*ElectricityProduced_MW*1000./1E6;
    
    else
        throw(MException('CapitalCost_SurfacePlant_ORC_GEOPHIRES:UnknownOrcModel','Unknown GEOPHIRES ORC Cost Model'));
    end
    
                
    %1.02 to convert cost from 2012 to 2016 #factor 1.15 for 15% contingency and 1.12 for 12% indirect costs.
    % Assume 2012 dollars
    % Geophires relations in millions of dollars
    C_plant_2012 = 1.12*1.15*Cplantcorrelation * 1e6; %*1.02;
    
    PPI_HX = PPI('PPI_HX', params.costYear) / PPI('PPI_HX', 2012);
    PPI_T_G = PPI('PPI_T-G', params.costYear) / PPI('PPI_T-G', 2012);
    PPI_PE = PPI('PPI_PE', params.costYear) / PPI('PPI_PE', 2012);
    PPI_avg = mean([PPI_HX, PPI_T_G, PPI_PE]);

    C_plant = C_plant_2012 * PPI_avg;
    
end

