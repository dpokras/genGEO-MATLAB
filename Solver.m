classdef Solver < handle
    properties
        func1 = [];
        allowExtrapolation = true;
        showPlot = false;
        creativeSolving = true;
        count = 0;
    end
    methods
        function obj = Solver(allowExtrapolation)
            if (nargin == 1)
                obj.allowExtrapolation = allowExtrapolation;
            end
        end
        % x is dependentVar, y is independentVar
        function x_zero = AddDataAndEstimate(obj, x, y)
            obj.count = obj.count + 1;
            obj.func1 = [obj.func1; x y];
            obj.func1 = sortrows(obj.func1);
            %remove duplicate points
            obj.func1 = unique(obj.func1,'rows');
            x_zero = Solver.ConvergeToZero(obj.func1, obj.allowExtrapolation, obj.showPlot);

            % if trouble converging, try to just guess half way between the
            % crossing point. Try once every 7 steps.
            if ( mod(obj.count,7)==0 && obj.creativeSolving == true)
                [row,~]=find(obj.func1(:,1)<x_zero);
                %it's a problem if nothing or all found
                if (isempty(row) || max(row)==size(row,1))
                    disp('Solver: no points on both sides of x_zero!');
                else
                    lastRowBelowXzero = max(row);
                    firstRowAboveXzero = max(row)+1;
                    % guess half way between
                    x1 = obj.func1(lastRowBelowXzero,1);
                    x2 = obj.func1(firstRowAboveXzero,1);
                    x_zero_new = (x1-x2)/2 + x2;
                    if (obj.showPlot == true)
                        disp(strcat(['Solver: Problems converging, old x_zero=' num2str(x_zero) ', new x_zero=' num2str(x_zero_new)]));    
                    end
                    x_zero = x_zero_new;
                end
            end
            
            % if debugging, show output
            if (obj.showPlot == true)
                disp(strcat(['Adding IV=' num2str(x) ' for DV=' num2str(y) ', predicting IV_0=' num2str(x_zero) '.']));
            end
        end
    end
    % Static Methods
    methods (Static)
        function x_zero = ConvergeToZero(func1, allowExtrapolation, showPlot)
            % This takes an array with 2 columns, x and y. It will guess the x value
            % needed to make y zero.

            % func2 is x-axis.
            func2 = [0, 0; 1, 0];

            try
                [x_zero, y_zero] = GetIntersection(func1, func2, allowExtrapolation);
            catch ME
                if (strcmp(ME.identifier,'GetIntersection:NoIntersection'))
                    % The lines appear divergent on either end.'
                    % Old-school nudge it in the right direction.
                    %x_zero = func1(end,1) - func1(end,2);
                    x_zero = NaN;
                    y_zero = 0;
                    disp(strcat(['Converge to zero not working. Nudging x']));
                elseif (strcmp(ME.identifier,'GetIntersection:NotEnoughRows'))
                    % Only single values
                    % if last y is zero, just return corresponding x.
                    if (func1(end,2)==0)
                        x_zero = func1(end,1);
                        y_zero = 0;
                    else
                        % Can't guess
                        x_zero = NaN;
                        y_zero = 0;
                    end
                else
                    rethrow(ME);
                end
            end

            %showPlot = true;
            if (showPlot == true)
                close(gcf);
                
                hold on;
                plot(func1(:,1),func1(:,2),'-o','DisplayName','func1');
                plot(func2(:,1),func2(:,2),'-o','DisplayName','func2');
                scatter(x_zero, y_zero, 100, 'pentagram', 'filled', 'm', 'DisplayName', 'Intersection');
                legend('FontSize',12);
                grid on;
                xlabel('independentVar');
                ylabel('dependentVar');
                hold off;
            end
        end
    end
end