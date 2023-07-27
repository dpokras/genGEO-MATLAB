function result = LCOE_Simple(CapitalCost, Capacity_W, Capacity_Q, params)

    % Capacity input in Watts
    % Capital cost input in $
    % LCOE output in $/W-h
    % SpCC output in $/W
    
    %Check heats & powers
    if (CapitalCost <= 0)
        result.LCOE = NaN;
        result.SpecificCapitalCost = NaN;
        return;
    elseif (params.F_OM < 0)
        throw(MException('LCOE_Simple:NegativeF_OM','Negative O&M Fraction')); 
    elseif (params.discountRate < 0)
        throw(MException('LCOE_Simple:NegativeDiscountRate','Negative Discount Rate')); 
    elseif (params.Lifetime <= 0)
        throw(MException('LCOE_Simple:NegativeLifetime','Negative Lifetime')); 
    elseif (params.CapacityFactor <= 0 || params.CapacityFactor > 1)
        throw(MException('LCOE_Simple:BadCapacityFactor','Bad Capacity Factor')); 
    end
    
    CRF = params.discountRate * (1+params.discountRate)^params.Lifetime / ((1+params.discountRate)^params.Lifetime - 1);
    
    result.LCOE = CapitalCost * (CRF + params.F_OM) / (Capacity_W * params.CapacityFactor * 8760);
    if Capacity_W > 0 && Capacity_Q == 0 
        result.SpCC_W = CapitalCost / Capacity_W;
        result.SpCC_Q = 0;
    elseif Capacity_Q > 0 && Capacity_W == 0 
        result.SpCC_Q = CapitalCost / Capacity_Q;
        result.SpCC_W = 0;
    else
        result.SpCC_Q = CapitalCost / Capacity_Q;
        result.SpCC_W = CapitalCost / Capacity_W;
    end
    result.SpCC_dH = CapitalCost / (Capacity_Q + Capacity_W);
end

