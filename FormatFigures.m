function [] = FormatFigures(fluid, Xtype, Ztype, legendLines, XLim, YLim, year, WACC, F_OM, dTdz, showProperties)

    %figurePosition = get(gcf,'Position');
    %set(gcf,'Position',[figurePosition(1),figurePosition(2),1000,900]);
    set(gcf,'Position',[50,50,1000,900]);
    set(gca,'FontSize',16);
    %set(gca,'XScale','log');
    set(gca,'XGrid','on');
    set(gca,'YGrid','on');
    if (strcmp(Xtype,'Transmissivity'))
        set(gca,'XLim',[2,XLim]);
        set(gca,'XTick',[2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]);
        set(gca,'XTickLabel',{'100','320','1 000','3 200','10 000','32 000','100 000','320 000','1 000 000','3 200 000','10 000 000'});
        xlabel('Reservoir Transmissivity [mD-m] (Log Scale)');
        % If XLim greater than 5, rotate axis labels so they fit
        if (XLim > 5)
            xtickangle(30);
        end
    elseif (strcmp(Xtype,'ResLength'))
        %set(gca,'XLim',[1000,10000]);
        set(gca,'XLim',[1000,XLim]);
        set(gca,'XTick',[1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000]);
        xlabel('Reservoir Horizontal Length [m]');
    end
    
    % Set Y axis limits
    set(gca,'YLim',[1000,YLim]);
    set(gca,'YTick',[1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000]);
    set(gca,'YTickLabel',[1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000]);
    ylabel('Reservoir Depth [m]');

    if (strcmp(Xtype,'Transmissivity'))
        [caption, titlestr] = FigureCaption(Ztype, fluid, legendLines, year, WACC, F_OM, dTdz);
    elseif (strcmp(Xtype,'ResLength'))
        [caption, titlestr] = FigureCaption_Conduction(Ztype, fluid, legendLines, year, WACC, F_OM);
    end
    title(titlestr);
    
    if (showProperties == true)
        annot = annotation('textbox','String',caption);
        %annot.Units = 'points';
        annot.Position = [0.15 0.6 0.3 0.3]; %[110 500 150 90];
        annot.FitBoxToText = 'on';
        annot.BackgroundColor = 'White';
        annot.FontSize = 11;
        %annot.Interpreter = 'latex';
    end
    
    legend('FontSize',12);
end

