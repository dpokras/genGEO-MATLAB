function [x_intercept, y_intercept] = GetIntersection(func1, func2, allowExtrapolation)

%Both functions have to have only two columns
if (size(func1, 2) ~= 2)
    throw(MException('GetIntersection:OnlyTwoColumns','Function1 does not have two columns'));
end
if (size(func2, 2) ~= 2)
    throw(MException('GetIntersection:OnlyTwoColumns','Function2 does not have two columns'));
end

%Both functions have to have at least two rows
if (size(func1, 1) < 2 || size(func2, 1) < 2)
    throw(MException('GetIntersection:NotEnoughRows','Functions do not both have more than two rows'));
end

% data points that directly cross
data_int = [];
% data points where the extrapolation of one cross another
extrap_data_int = [];
% data points where two extrapolations cross
extrap_int = [];

% Iterate through all possible line combinations
for index1=1:(size(func1,1)-1)
    for index2=1:(size(func2,1)-1)
        % Now x2(index2) is between x1(index1 -1) and x1(index1)
        x1a = func1(index1,1);
        y1a = func1(index1,2);
        x1b = func1(index1+1,1);
        y1b = func1(index1+1,2);
        x2a = func2(index2,1);
        y2a = func2(index2,2);
        x2b = func2(index2+1,1);
        y2b = func2(index2+1,2);

        % Solve for the equation of each line
        m1 = (y1b - y1a) / (x1b - x1a);
        b1 = y1a - m1*x1a;
        m2 = (y2b - y2a) / (x2b - x2a);
        b2 = y2a - m2*x2a;
        % Solve for x where the lines intersect
        x_int = (b2 - b1) / (m1 - m2);
        y_int = m1*x_int + b1;
        % If x_int is between x1a and x1b, it crosses here
        if (x_int >= x1a && x_int < x1b && x_int >= x2a && x_int < x2b)
            % Found an intersection!
            data_int = [data_int; x_int y_int];
        elseif (allowExtrapolation == true)
            % See if extrapolated ends intersect anywhere
            if (index2 == 1 && sign(x_int-x2a)==sign(x2a-x2b) && x_int >= x1a && x_int < x1b)
                %It intersects this extrapolation!
                extrap_data_int = [extrap_data_int; x_int y_int];
            elseif (index2 == size(func2,1)-1 && sign(x_int-x2b)==sign(x2b-x2a) && x_int >= x1a && x_int < x1b)
                %It intersects the extrapolation!
                extrap_data_int = [extrap_data_int; x_int y_int];
            end
            
            if (index1 == 1)
                if (sign(x_int-x1a)==sign(x1a-x1b) || size(func1,1)==2)
                    if (x_int >= x2a && x_int < x2b)
                        %it intersects the existing line!
                        extrap_data_int = [extrap_data_int; x_int y_int];
                    elseif (index2 == 1 && (sign(x_int-x2a)==sign(x2a-x2b) || size(func2,1)==2))
                        %It intersects this extrapolation!
                        extrap_int = [extrap_int; x_int y_int];
                    elseif (index2 == size(func2,1)-1 && (sign(x_int-x2b)==sign(x2b-x2a) || size(func2,1)==2))
                        %It intersects the extrapolation!
                        extrap_int = [extrap_int; x_int y_int];
                    end
                end
            elseif (index1 == size(func1,1)-1)
                if (sign(x_int-x1b)==sign(x1b-x1a))
                    if (x_int >= x2a && x_int < x2b)
                        %it intersects the existing line!
                        extrap_data_int = [extrap_data_int; x_int y_int];
                    elseif (index2 == 1 && (sign(x_int-x2a)==sign(x2a-x2b) || size(func2,1)==2))
                        %It intersects this extrapolation!
                        extrap_int = [extrap_int; x_int y_int];
                    elseif (index2 == size(func2,1)-1 && (sign(x_int-x2b)==sign(x2b-x2a) || size(func2,1)==2))
                        %It intersects the extrapolation!
                        extrap_int = [extrap_int; x_int y_int];
                    end
                end
            end
        end
    end
end

% Now we choose which of the results we want
if (~isempty(data_int))
    % take first
    x_intercept = data_int(1,1);
    y_intercept = data_int(1,2);
elseif (~isempty(extrap_data_int))
    % take first
    x_intercept = extrap_data_int(1,1);
    y_intercept = extrap_data_int(1,2);
elseif (~isempty(extrap_int))
    % find closest value
    x_dx_min = Inf;
    for r=1:size(extrap_int,1)
        x_int_left_dx = max([func1(:,1); func2(:,1)]) - extrap_int(r,1);
        if (x_int_left_dx < 0)
            x_int_left_dx = NaN;
        end
        x_int_right_dx = extrap_int(r,1) - min([func1(:,1); func2(:,1)]);
        if (x_int_right_dx < 0)
            x_int_right_dx = NaN;
        end
        x_int_dx = min(x_int_left_dx,x_int_right_dx);
        if (x_int_dx < x_dx_min)
            x_dx_min = x_int_dx;
            % use this one
            x_intercept = extrap_int(r,1);
            y_intercept = extrap_int(r,2);
        end
    end
    %make sure x_dx_min was set
    if (isinf(x_dx_min))
        throw(MException('GetIntersection:NoIntersection','Functions do not intersect'));
    end
else
    % Couldn't find any intersections
    throw(MException('GetIntersection:NoIntersection','Functions do not intersect'));
end
    

end

