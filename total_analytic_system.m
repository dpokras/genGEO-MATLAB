function result = total_analytic_system(params)

    if (strcmp(params.fluid,'Water'))
        if (strcmp(params.simulatorType,'genGEO'))
            result = total_analytic_system_water(params);
        elseif (strcmp(params.simulatorType,'GEOPHIRES'))
            result = total_analytic_system_GEOPHIRES(params);
        else
            throw(MException('totalSystem:NotImplemented','simulatorType Not Implemented'));
        end
    elseif (strcmp(params.fluid,'CO2'))
        result = total_analytic_system_co2(params);
    else
        throw(MException('totalSystem:NotImplemented','Fluid Not Implemented'));
    end
    
end

