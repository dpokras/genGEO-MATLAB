clear;
close all;
clc;

%% Define System Parameters
params = SimulationParameters; %radius, surface temp, dT/dz, etc.
params.wellFieldType = 'Doublet'; 
params.time_years = 30;
params.dT_dz = 0.035;
params.optimizationMode = 'MinimizeLCOE_Greenfield';
params.orcFluid = 'R245fa';
params.well_radius = 0.25;

for b = 1:2
    if b == 1
        base_case = "Eavor"; % 2 Power Plants each with 4RW, 1IW, 1PW
    elseif b == 2
        base_case = "Deep Well";
    end
%% Define parameters for base cases  
    if base_case == "Eavor"
        params.res_length = 5000; % in meters
        params.depth = 3500; % in meters
        depth1 = params.depth;
        params.system = 'Conduction4';
    elseif base_case == "Deep Well"
        params.res_length = 5000; % in meters
        params.depth = 8000; % in meters
        depth2 = params.depth;
        params.system = 'Conduction4';
    end
    fluid_All = ["Water","CO2"];
    %% For Each Fluid, first find optimum m-dot, then find temperature and pressure profiles
    for j = 1:2
        params.fluid = fluid_All(j);
        % Determine Optimum m-dot for fluid
        %result = total_analytic_system_optmdot(params);
        params.m_dot_IP = 22;
        time_seconds = 3600 * 24 * 365 * params.time_years; %seconds
        if b == 1
            temp_at_bottom = params.T_surface_rock_C + params.dT_dz * depth1; %C
        elseif b == 2
            temp_at_bottom = params.T_surface_rock_C + params.dT_dz * depth2; %C
        end
        if params.fluid == "Water"
            result = total_analytic_system_water(params);
            z_IW_Water = result.injection.State(:,1)*(-1);
            z_LW_Water = (result.reservoir.State(:,1)- params.depth)*(-1);
            z_PWl_Water = (result.productionLower.State(:,1) - params.depth)*(-1);
            z_PWu_Water = result.productionUpper.State(:,1);
            z_PWu_Water = flip(z_PWu_Water);
            T_IW_Water = result.injection.State(:,3);
            T_LW_Water = result.reservoir.State(:,3);
            T_PWl_Water = result.productionLower.State(:,3);
            T_PWu_Water = result.productionUpper.State(:,3);
            P_IW_Water = result.injection.State(:,2);
            P_LW_Water = result.reservoir.State(:,2);
            P_PWl_Water = result.productionLower.State(:,2);
            P_PWu_Water = result.productionUpper.State(:,2);
        end
        
        
        if params.fluid == "CO2"
            result = total_analytic_system_co2(params);
            z_IW_CO2 = result.injection.State(:,1)*(-1);
            z_LW_CO2 = (result.reservoir.State(:,1)- params.depth)*(-1);
            z_PWu_CO2 = result.productionUpper.State(:,1);
            z_PWu_CO2 = flip(z_PWu_CO2);
            T_IW_CO2 = result.injection.State(:,3);
            T_LW_CO2 = result.reservoir.State(:,3);
            T_PWu_CO2 = result.productionUpper.State(:,3);
            P_IW_CO2 = result.injection.State(:,2);
            P_LW_CO2 = result.reservoir.State(:,2);
            P_PWu_CO2 = result.productionUpper.State(:,2);
        end
        dT_dz_plot = params.T_surface_rock_C + z_IW_Water.* params.dT_dz;
        %% Plot Figure
        figure(1)
        if b == 1
            subplot(1,2,1)
                dTdz1 = plot(dT_dz_plot,z_IW_Water, '-','Color',"#a8a8a8",'linewidth',1.4);
                hold on
                if params.fluid == "Water" && b == 1
                    w_T_inj = plot(T_IW_Water,z_IW_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                    w_P_inj = plot(P_IW_Water,z_IW_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                    w_T_res = plot(T_LW_Water,z_LW_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                    w_P_res = plot(P_LW_Water,z_LW_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                    w_T_prodl = plot(T_PWl_Water,z_PWl_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                    w_T_produ = plot(T_PWu_Water,z_PWu_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                    w_P_prodl = plot(P_PWl_Water,z_PWl_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                    w_P_produ = plot(P_PWu_Water,z_PWu_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                    %Plot labels for T-profile
                    plot(T_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(1)));
                    text(T_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2' );
                    plot(T_LW_Water(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(2)));
                    text(T_LW_Water(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                    plot(T_LW_Water(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(3)));
                    text(T_LW_Water(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                    plot(T_PWu_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(4)));
                    text(T_PWu_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2');
                    %Plot labels for P-profile
                    plot(P_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(1)));
                    text(P_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2' );
                    plot(P_LW_Water(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(2)));
                    text(P_LW_Water(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                    plot(P_LW_Water(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(3)));
                    text(P_LW_Water(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                    plot(P_PWu_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(4)));
                    text(P_PWu_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2');
                    sIW_water = [z_IW_Water,T_IW_Water,P_IW_Water];
                    sLW_water = [z_LW_Water,T_LW_Water,P_LW_Water];
                    sPWl_water = [z_PWl_Water,T_PWl_Water,P_PWl_Water];
                    sPWu_water = [z_PWu_Water,T_PWu_Water,P_PWu_Water];
                    table_water_shallow = [sIW_water,sLW_water,sPWl_water,sPWu_water];
                    table_water_shallow = array2table(table_water_shallow)
                end
                if params.fluid == "CO2" && b == 1
                    c_T_inj = plot(T_IW_CO2,z_IW_CO2,'-k','linewidth',1.4);
                    c_P_inj = plot(P_IW_CO2,z_IW_CO2,'--k','linewidth',1.4);
                    c_T_res = plot(T_LW_CO2,z_LW_CO2,'-k','linewidth',1.4);
                    c_P_res = plot(P_LW_CO2,z_LW_CO2,'--k','linewidth',1.4);
                    c_T_produ = plot(T_PWu_CO2,z_PWu_CO2,'-k','linewidth',1.4);
                    c_P_produ = plot(P_PWu_CO2,z_PWu_CO2,'--k','linewidth',1.4);
                    %Plot labels for T-profile
                    plot(T_IW_CO2(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(1)));
                    text(T_IW_CO2(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k' );
                    plot(T_LW_CO2(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(2)));
                    text(T_LW_CO2(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                    plot(T_LW_CO2(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(3)));
                    text(T_LW_CO2(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                    plot(T_PWu_CO2(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(4)));
                    text(T_PWu_CO2(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k');
                    %Plot labels for P-profile
                    plot(P_IW_CO2(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(1)));
                    text(P_IW_CO2(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k' );
                    plot(P_LW_CO2(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(2)));
                    text(P_LW_CO2(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                    plot(P_LW_CO2(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(3)));
                    text(P_LW_CO2(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                    plot(P_PWu_CO2(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                    strValues = strtrim(cellstr(num2str(4)));
                    text(P_PWu_CO2(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k');
                    sIW_co2 = [z_IW_CO2,T_IW_CO2,P_IW_CO2];
                    sLW_co2 = [z_LW_CO2,T_LW_CO2,P_LW_CO2];
                    sPWu_co2 = [z_PWu_CO2,T_PWu_CO2,P_PWu_CO2];
                    table_co2_shallow = [sIW_co2,sLW_co2,sPWu_co2];
                    table_co2_shallow = array2table(table_co2_shallow);
                end
                set(gca, 'YDir','reverse')
                set(gca,'FontSize',11)
                xlim([0,temp_at_bottom]) % Specifying limits on x-axis
                xlabel({'Temperature [\circ{}C]', 'Pressure [MPa]'}, 'FontSize',11)  % x-axis label
                ylabel('Vertical Well Depth [m]', 'FontSize',11) % y-axis label
                title({'Temperature, Pressure Profiles', 'for 3.5 km Vertical Well Depth'}, 'FontSize',11)
                if params.fluid == "CO2" && b == 1
                    lgd = legend([dTdz1 w_T_inj w_P_inj c_T_inj c_P_inj],'Geothermal Gradient','Water Temperature','Water Pressure', 'CO_2 Temperature','CO_2 Pressure');
                    lgd.Title.String = 'Legend';
                    lgd.Title.FontSize = 11;
                    lgd.FontSize = 11;
                end
                grid on

%         resevoir_axis = linspace(0,res_length,101);
     elseif b == 2
        subplot(1,2,2)
            dTdz2 = plot(dT_dz_plot,z_IW_Water, '-','Color',"#a8a8a8",'linewidth',1.4);
            hold on
            if params.fluid == "Water" && b == 2
                w_T_inj = plot(T_IW_Water,z_IW_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                w_P_inj = plot(P_IW_Water,z_IW_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                w_T_res = plot(T_LW_Water,z_LW_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                w_P_res = plot(P_LW_Water,z_LW_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                w_T_prodl = plot(T_PWl_Water,z_PWl_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                w_T_produ = plot(T_PWu_Water,z_PWu_Water,'-','Color',"#4e5fa2",'linewidth',1.4);
                w_P_prodl = plot(P_PWl_Water,z_PWl_Water,'--','Color',"#4e5fa2",'linewidth',1.4);
                w_P_produ = plot(P_PWu_Water,z_PWu_Water,'--','Color',"#4e5fa2",'linewidth',1.4);

                %Plot labels for T-profile
                plot(T_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(1)));
                text(T_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2' );
                plot(T_LW_Water(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(2)));
                text(T_LW_Water(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                plot(T_LW_Water(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(3)));
                text(T_LW_Water(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                plot(T_PWu_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(4)));
                text(T_PWu_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2');
                %Plot labels for P-profile
                plot(P_IW_Water(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(1)));
                text(P_IW_Water(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2' );
                plot(P_LW_Water(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(2)));
                text(P_LW_Water(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                plot(P_LW_Water(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(3)));
                text(P_LW_Water(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', '#4e5fa2');
                plot(P_PWu_Water(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','#4e5fa2', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(4)));
                text(P_PWu_Water(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', '#4e5fa2');
                dIW_water = [z_IW_Water,T_IW_Water,P_IW_Water];
                dLW_water = [z_LW_Water,T_LW_Water,P_LW_Water];
                dPWl_water = [z_PWl_Water,T_PWl_Water,P_PWl_Water];
                dPWu_water = [z_PWu_Water,T_PWu_Water,P_PWu_Water];
                table_water_deep = [dIW_water,dLW_water,dPWl_water,dPWu_water];
                table_water_deep = array2table(table_water_deep)
            end
            if params.fluid == "CO2" && b == 2
                c_T_inj = plot(T_IW_CO2,z_IW_CO2,'-k','linewidth',1.4);
                c_P_inj = plot(P_IW_CO2,z_IW_CO2,'--k','linewidth',1.4);
                c_T_res = plot(T_LW_CO2,z_LW_CO2,'-k','linewidth',1.4);
                c_P_res = plot(P_LW_CO2,z_LW_CO2,'--k','linewidth',1.4);
                c_T_produ = plot(T_PWu_CO2,z_PWu_CO2,'-k','linewidth',1.4);
                c_P_produ = plot(P_PWu_CO2,z_PWu_CO2,'--k','linewidth',1.4);
                %Plot labels for T-profile
                plot(T_IW_CO2(1),0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(1)));
                text(T_IW_CO2(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k' );
                plot(T_LW_CO2(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(2)));
                text(T_LW_CO2(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                plot(T_LW_CO2(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(3)));
                text(T_LW_CO2(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                plot(T_PWu_CO2(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(4)));
                text(T_PWu_CO2(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k');
                %Plot labels for P-profile
                plot(P_IW_CO2(1), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(1)));
                text(P_IW_CO2(1),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k' );
                plot(P_LW_CO2(1), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(2)));
                text(P_LW_CO2(1),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                plot(P_LW_CO2(end), params.depth, 'ok', 'MarkerSize', 10, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(3)));
                text(P_LW_CO2(end),params.depth,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8,'Color', 'k');
                plot(P_PWu_CO2(end), 0, 'ok', 'MarkerSize', 10,'MarkerEdgeColor','k', 'MarkerFaceColor', 'w'); grid on; hold on
                strValues = strtrim(cellstr(num2str(4)));
                text(P_PWu_CO2(end),0,strValues, 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'fontsize', 8, 'Color', 'k');
                dIW_co2 = [z_IW_CO2,T_IW_CO2,P_IW_CO2];
                dLW_co2 = [z_LW_CO2,T_LW_CO2,P_LW_CO2];
                dPWu_co2 = [z_PWu_CO2,T_PWu_CO2,P_PWu_CO2];
                table_co2_deep = [dIW_co2,dLW_co2,dPWu_co2];
                table_co2_deep = array2table(table_co2_deep)
            end
            set(gca, 'YDir','reverse')
            set(gca,'FontSize',11)
            xlim([0,temp_at_bottom]) % Specifying limits on x-axis
            xlabel({'Temperature [\circ{}C]', 'Pressure [MPa]'}, 'FontSize',11)  % x-axis label
            ylabel('Vertical Well Depth [m]', 'FontSize',11) % y-axis label
            title({'Temperature, Pressure Profiles', 'for 8.0 km Vertical Well Depth'}, 'FontSize',11)
            grid on    
        end
    end

%     if params.fluid == "CO2" && b == 1
%         sIW_co2 = [z_IW_CO2;T_IW_CO2;P_IW_CO2];
%         sLW_co2 = [z_LW_CO2;T_LW_CO2;P_LW_CO2];
%         sPWu_co2 = [z_PWu_CO2;T_PWu_CO2;P_PWu_CO2];
%         table_co2_shallow = [sIW_co2;sLW_co2;sPWu_co2];
%         table_co2_shallow = array2table(table_co2_shallow)
%     end
%     if params.fluid == "Water" && b == 2
%         dIW_water = [z_IW_Water,T_IW_Water,P_IW_Water];
%         dLW_water = [z_LW_Water,T_LW_Water,P_LW_Water];
%         dPWl_water = [z_PWl_Water,T_PWl_Water,P_PWl_Water];
%         dPWu_water = [z_PWu_Water,T_PWu_Water,P_PWu_Water];
%         table_water_deep = [dIW_water,dLW_water,dPWl_water,dPWu_water];
%         table_water_deep = array2table(table_water_deep)
%     end
%     if params.fluid == "CO2" && b == 2
%         dIW_co2 = [z_IW_CO2,T_IW_CO2,P_IW_CO2];
%         dLW_co2 = [z_LW_CO2,T_LW_CO2,P_LW_CO2];
%         dPWu_co2 = [z_PWu_CO2,T_PWu_CO2,P_PWu_CO2];
%         table_co2_deep = [dIW_co2,dLW_co2,dPWu_co2];
%         table_co2_deep = array2table(table_co2_deep)
%     end
end

writetable(table_water_shallow,'data\Figure9_Data_TP_Profiles_water_3_5.xlsx')
writetable(table_co2_shallow,'data\Figure9_Data_TP_Profiles_co2_3_5.xlsx')
writetable(table_water_deep,'data\Figure9_Data_TP_Profiles_water_8.xlsx')
writetable(table_co2_deep,'data\Figure9_Data_TP_Profiles_co2_8.xlsx')

caption = {'{\bfProperties}'...
    ,'Massflow to Minimize Cost'...
    ,'5000 m Lateral Well Length'...
    ,'35 �C km^{-1} Geothermal Gradient'...
    ,'15 �C Surface Temp'...
    ,'0.5 m Well Diameter'...
    ,'Four Horizontal Laterals'...
    ,strcat(['Year-' num2str(30) ' Depletion Values'])...
    ,'Malek|Adams|Rossi|Schiegg|Saar (2021)'};
annot = annotation('textbox','String',caption);
annot.Position = [0.23 0.6 0.25 0.3];
annot.FitBoxToText = 'on';
annot.BackgroundColor = 'White';
annot.FontSize = 11;
        
annota = annotation('textbox',[0.08 0.68 0.05 0.05],'String',"(a)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');
annotb = annotation('textbox',[0.55 0.68 0.05 0.05],'String',"(b)",'EdgeColor','none','FontSize',14,'fontweight', 'bold');

saveas(figure(1),'images\T_P_Profiles.fig');
