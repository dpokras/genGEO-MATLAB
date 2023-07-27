function [X,Y,Z,Z_Combined_Mask] = GetContourData(fluid, useBest, Xtype, Ztype, xlsxFile, co2Sheet, waterSheet)
    % Depth column 2
    % ResLength column 3
    % Transmissivity column 24
    % Power column 25
    % SpCC column 26
    % LCOE column 27
    depthColumn = 2;
    resLengthColumn = 3;
    transColumn = 46;
    powerColumn = 47;
    spccBrownfieldColumn = 48;
    spccGreenfieldColumn = 49;
    lcoeBrownfieldColumn = 50;
    lcoeGreenfieldColumn = 51;
    
    % Xtype
    if (strcmp(Xtype,'Transmissivity'))
        Xcolumn = transColumn;
    elseif (strcmp(Xtype,'ResLength'))
        Xcolumn = resLengthColumn;
    else
        throw(MException('GetContourData:XTypeNotImplemented','XType Not Implemented'));
    end
    
    % Ytype
    Ycolumn = depthColumn;
    
    % Ztype
    if (strcmp(Ztype,'Power'))
        Zcolumn = powerColumn;
    elseif (strcmp(Ztype,'SpCC_Brownfield'))
        Zcolumn = spccBrownfieldColumn;
    elseif (strcmp(Ztype,'SpCC_Greenfield'))
        Zcolumn = spccGreenfieldColumn;
    elseif (strcmp(Ztype,'LCOE_Brownfield'))
        Zcolumn = lcoeBrownfieldColumn;
    elseif (strcmp(Ztype,'LCOE_Greenfield'))
        Zcolumn = lcoeGreenfieldColumn;
    else
        throw(MException('GetContourData:ZTypeNotImplemented','ZType Not Implemented'));
    end
    

    if (useBest)
        %Use Both Datasets
        a_CO2 = readmatrix(xlsxFile,'Sheet',co2Sheet);
        b_Water = readmatrix(xlsxFile,'Sheet',waterSheet);
        sizeA = size(a_CO2,1);
        sizeB = size(b_Water,1);
        %Depth
        y = unique([a_CO2(:,Ycolumn); b_Water(:,Ycolumn)]);
        %Transmissivity
        x = unique([a_CO2(:,Xcolumn); b_Water(:,Xcolumn)]);
    elseif (strcmp(fluid,'CO2'))
        a_CO2 = readmatrix(xlsxFile,'Sheet',co2Sheet);
        sizeA = size(a_CO2,1);
        y = unique(a_CO2(:,Ycolumn));
        x = unique(a_CO2(:,Xcolumn));
    elseif (strcmp(fluid,'Water'))
        b_Water = readmatrix(xlsxFile,'Sheet',waterSheet);
        sizeB = size(b_Water,1);
        y = unique(b_Water(:,Ycolumn));
        x= unique(b_Water(:,Xcolumn));
    else
        throw(MException('getContourData:NotImplemented','Not Implemented'));
    end

    
    if (strcmp(Xtype,'Transmissivity'))
        % Use a log scale
        x = log10(x);
    elseif (strcmp(Xtype,'ResLength'))
        % Use normal scale
    end
    
    % Remove Missing NaN
    x = rmmissing(x);
    y = rmmissing(y);
    
    [X,Y] = meshgrid(x,y);

    sizeX = size(x,1);
    sizeY = size(y,1);
    Z_CO2 = NaN(sizeY,sizeX);
    Z_Water = NaN(sizeY,sizeX);
    Z_Combined = NaN(sizeY,sizeX);
    Z_Combined_Mask = NaN(sizeY,sizeX);

    %CO2
    if (useBest || strcmp(fluid,'CO2'))
        for i = 1:sizeA
            depth = a_CO2(i,Ycolumn);
            trans = a_CO2(i,Xcolumn);
            power = a_CO2(i,Zcolumn);

            [depth_tf,depth_index] = ismember(depth,y);
            if (strcmp(Xtype,'Transmissivity'))
                % log scale
                [trans_tf,trans_index] = ismember(log10(trans),x);
            elseif (strcmp(Xtype,'ResLength'))
                % Use normal scale
                [trans_tf,trans_index] = ismember(trans,x);
            end
            
            if (depth_tf == 0 || trans_tf == 0)
                continue;
            end

            Z_CO2(depth_index,trans_index) = power;
        end
    end
    %Water
    if (useBest || strcmp(fluid,'Water'))
        for i = 1:sizeB
            depth = b_Water(i,Ycolumn);
            trans = b_Water(i,Xcolumn);
            power = b_Water(i,Zcolumn);

            [depth_tf,depth_index] = ismember(depth,y);
            if (strcmp(Xtype,'Transmissivity'))
                [trans_tf,trans_index] = ismember(log10(trans),x);
            elseif (strcmp(Xtype,'ResLength'))
                [trans_tf,trans_index] = ismember(trans,x);
            end


            if (depth_tf == 0 || trans_tf == 0)
                continue;
            end

            Z_Water(depth_index,trans_index) = power;
        end
    end

    % Make some NaN if useBest
    if (useBest)
        for i = 1:sizeY
            for j = 1:sizeX
                if (isnan(Z_CO2(i,j)) && ~isnan(Z_Water(i,j)))
                    Z_Combined(i,j) = Z_Water(i,j);
                    continue;
                end
                if (isnan(Z_Water(i,j)) && ~isnan(Z_CO2(i,j)))
                    Z_Combined(i,j) = Z_CO2(i,j);
                    Z_Combined_Mask(i,j) = 1;
                    continue;
                end
                if (Z_CO2(i,j) > Z_Water(i,j))
                    Z_Water(i,j) = NaN;
                    Z_Combined(i,j) = Z_CO2(i,j);
                    Z_Combined_Mask(i,j) = 1;
                else
                    Z_CO2(i,j) = NaN;
                    Z_Combined(i,j) = Z_Water(i,j);
                end
            end
        end
    end

    if (useBest)
        Z = Z_Combined;
    elseif (strcmp(fluid,'CO2'))
        Z = Z_CO2;
    elseif (strcmp(fluid,'Water'))
        Z = Z_Water;
    end
end

