
assignment_(@eulerBackward)

function assignment_c(e_method)
    f = @(t, y) sin(3*t)- 2*y;  % Anonymus function for dy/dt
    t = 0; T = 8;               % Start and end values of interval
    n = [50 100 200 400];       % The diffrent values of n to determine
    y = 6/5;                    % Initial value of function 

    y_end = (93/65)*exp(-2*T)+(2/13)*sin(3*T)-(3/13)*cos(3*T);

    errors = zeros(length(n),2);
    i = 1;
    for ni=n
        a = e_method(t,T,ni,f,y);       % Computes points of equation
        errors(i, :) = [(T-t)/ni abs(y_end - a(end,2))];
        i = i + 1; 
    end

    loglog(errors(:, 1), errors(:,2), 'o')
end

% Function to get results from assignment 2.b
function assignment_b(e_method)
    f = @(t, y) sin(3*t)- 2*y;  % Anonymus function for dy/dt
    t = 0; T = 8;               % Start and end values of interval
    n = [50 100 200 400];       % The diffrent values of n to determine
    y = 6/5;                    % Initial value of function 
    
    % For every resolution n, compute the points of the differential
    % equation and plot the points. 
    for ni=n
        a = e_method(t,T,ni,f,y)        % Computes points of equation
        hold on
        plot(a(:,1),a(:,2), '--')        % Plots points
        hold off
    end
    
    hold on
    % Plot the exact solution
    x_plot = t:0.01:T;
    exact_solution = (93/65)*exp(-2*x_plot)+(2/13)*sin(3*x_plot)-(3/13)*cos(3*x_plot);
    plot(x_plot,exact_solution)
    
    % Add lables to functions
    legend('n = 50', 'n = 100', 'n = 200', 'n = 400', 'exact solution')
    hold off
end


% Function to compute approximate points of differential equation
%   params: t (starting value of interval), T (end value of interval),
%       n (number of discrete steps), f (anonymus function where dy/dt=f(t,y)), 
%       y (initial value y0)
function [out] = eulerForward(t, T, n, f, y)
    h = (T-t)/n;

    out = zeros(n+1, 2);
    out(1, :) = [t y];
    for i=1:n
        y = y + h*f(t, y);
        t = t + h;
        out(i+1, :) = [t y];
    end
end

function [out] = eulerBackward(t, T, n, f, y)
    h = (T-t)/n;

    f = @(ta, ya) (ya+h*sin(3*ta))/(1+2*h)

    out = zeros(n+1, 2);
    out(1, :) = [t y];
    for i=1:n
        y = f(t, y);
        t = t + h;
        out(i+1, :) = [t y];
    end
end
