classdef CoolProp
    methods(Static)
        function outPropVal = PropsSI(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid)
        % This function makes a coolprop call. In this case of CO2 at pressures
        % very near the critical point, it will do some linear interpolation from
        % nearby values.
        
            if numel(size(prop1Val)) >= 3 || numel(size(prop2Val)) >= 3  
                
                %flatten array to fit into CoolProp
                prop1Val_flat = reshape(prop1Val, 1, prod(size(prop1Val)));
                prop2Val_flat = reshape(prop2Val, 1, prod(size(prop2Val)));


                prop1Val_cell = mat2cell(prop1Val_flat,1,ones(1,size(prop1Val_flat,2)));
                prop2Val_cell = mat2cell(prop2Val_flat,1,ones(1,size(prop2Val_flat,2)));
                outPropVal_cell = cell(py.CoolProp.CoolProp.PropsSI(outPropName, prop1Name, prop1Val_cell, prop2Name, prop2Val_cell, fluid));
                outPropVal = cellfun(@double, outPropVal_cell);

                outPropVal = reshape(outPropVal, size(prop1Val,1),size(prop1Val,2),size(prop1Val,3));
            
            elseif numel(size(prop1Val)) >= 2 || numel(size(prop2Val)) >= 2
                
                prop1Val_cell = mat2cell(prop1Val,1,ones(1,size(prop1Val,2)));
                prop2Val_cell = mat2cell(prop2Val,1,ones(1,size(prop2Val,2)));
                outPropVal_cell = cell(py.CoolProp.CoolProp.PropsSI(outPropName, prop1Name, prop1Val_cell, prop2Name, prop2Val_cell, fluid));
                outPropVal = cellfun(@double, outPropVal_cell);
            else
                outPropVal = py.CoolProp.CoolProp.PropsSI(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid);
            end
        end

        function outPhase = PhaseSI(prop1Name, prop1Val, prop2Name, prop2Val, fluid)
        % This function makes a coolprop call to identify the phase of the
        % fluid at that state.

            outPhase = py.CoolProp.CoolProp.PhaseSI(prop1Name, prop1Val, prop2Name, prop2Val, fluid);
        end
    end
end
% if (strcmp(prop1Name,'P'))
%     % pressure is first property
%     P = prop1Val;
%     xN = prop2Name;
%     xV = prop2Val;
% elseif (strcmp(prop2Name,'P'))
%     % pressure is second property
%     P = prop2Val;
%     xN = prop1Name;
%     xV = prop1Val;
% else
%     % no pressure at all
%     outPropVal = CoolProp_Library(outPropName, prop1Name, prop1Val, prop2Name, prop2Val, fluid);
%     return;
% end
% 
% %P_crit_co2 = CoolProp('PCRIT', "", 0, "", 0, 'CO2');
% P_crit_co2 = 7377300; %7.3773e6;
% dP_tolerance = 5e3; % 500 kPa
% P_upper = P_crit_co2 + dP_tolerance;
% P_lower = P_crit_co2 - dP_tolerance;
% 
% if (strcmp(fluid,'CO2') && P > P_lower && P < P_upper )
%     % If pressure is very near and below the critical point of 7.377 MPa (for CO2), there will be a coolprop convergence error 
%     %disp(strcat(['Semi_analytic_well:Manually adjusting pressure from ' num2str(P/1e6) ' MPa to 7.37 MPa to avoid CoolProp CO2 critical point convergence issues.']));
%     % linear interpolate between two extremes
%     outPropVal_lower = CoolProp_Library(outPropName, 'P', P_lower, xN, xV, fluid);
%     outPropVal_upper = CoolProp_Library(outPropName, 'P', P_upper, xN, xV, fluid);
%     outPropVal = interp1([P_lower, P_upper], [outPropVal_lower, outPropVal_upper], P);
% else
%     outPropVal = CoolProp_Library(outPropName, 'P', P, xN, xV, fluid);
% end
% 
% end
% 
% 
% % this function makes sure the appropriate coolprop library is called. On
% % windows machines this calls the python library. On macs, I wasn't able to
% % get the python library to work, so it calls an older coolprop matlab
% % wrapper.
% function [outValue] = CoolProp_Library(outName, in1Name, in1Value, in2Name, in2Value, fluid)
% %addpath('/Users/badams/polybox/MacExchange/Matlab_Coolprop/Coolprop');
% %addpath('/Users/badams/polybox/MacExchange/Matlab_Coolprop');
% %CoolProp http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
% 
% 
% try
%     outValue = py.CoolProp.CoolProp.PropsSI(outName, in1Name, in1Value, in2Name, in2Value, fluid);
% catch ME
%     if (strcmp(ME.identifier,'MATLAB:undefinedVarOrClass'))
%         %disp(strcat(['CoolProp:Didnt find coolprop.']));
% 
%         %path = '/Users/badams/Matlab_Coolprop';
%         %userPath = userpath;
%         %onPath = contains(lower(userPath),lower(path));
%         %if (onPath == false)
%             %addpath('/Users/badams/Matlab_Coolprop/Coolprop');
%             %addpath('/Users/badams/Matlab_Coolprop');
%             %onPath = ~isempty(strfind(lower(userpath),lower(path)));
%         %end
%         %%outValue = PropsSI('T','P',100000,'Q',0,'Water');
%         outValue = PropsSI(outName, in1Name, in1Value, in2Name, in2Value, fluid);
%     else
%         rethrow(ME)
%     end
% end
% end

