                                %% Project phase 1 %%

function main()

    range = input('Enter the range [a, b]: ');
    fprintf('Searching in the range (%.3f, %.3f)...\n', range(1), range(2));
    
    Exhaustive_Search(range);
    Bisection_Method(range);
end

function Exhaustive_Search(range)
    a = range(1);
    b = range(2);
    n = input('Enter the number of steps (n): ');
    delta = (b - a) / n;
    x1 = a;
    x2 = x1 + delta;
    x3 = x2 + delta;
    fx1 = objective_function(x1);
    fx2 = objective_function(x2);
    fx3 = objective_function(x3);
    funeval = 3;
    
    fprintf('\n==============Exhaustive Search=============\n');
    fprintf('#It\t\tx1\t\tx2\t\tx3\t\tf(x1)\t\tf(x2)\t\tf(x3)\n');
    
    while x3 <= b
        fprintf('%d\t%.3f\t%.3f\t%.3f\t%10.3f\t%10.3f\t%f\n', funeval, x1, x2, x3, fx1, fx2, fx3);
        
        if (fx1 >= fx2 && fx2 <= fx3) || (fx1 <= fx2 && fx2 >= fx3)
            break;
        end
        
        x1 = x2;
        x2 = x3;
        x3 = x2 + delta;
        fx1 = fx2;
        fx2 = fx3;
        fx3 = objective_function(x3);
        funeval = funeval + 1;
    end
    
    extrema = 'minimum';
    if fx1 > fx2
        extrema = 'maximum';
    end
    
    fprintf('The %s point lies between (%.3f, %.3f)\n', extrema, x1, x3);
    fprintf('Total number of function evaluations: %d\n', funeval);
    fprintf('============================================\n\n');
    range(1) = x1;
    range(2) = x3;
end

function Bisection_Method(range)
    a = range(1);
    b = range(2);
    x1 = a;
    x2 = b;
    dfx1 = first_derivative_function(x1);
    dfx2 = first_derivative_function(x2);
    epsilon = 0.001;
    funeval = 0;
    
    fprintf('\n==============Bisection Method==============\n');
    
    if dfx1 * dfx2 >= 0 && dfx1 - dfx2 > epsilon
        fprintf('The input range is not appropriate for the bisection method\n');
        return;
    elseif dfx1 > 0 && dfx2 < 0
        x1 = b;
        x2 = a;
        temp = dfx1;
        dfx1 = dfx2;
        dfx2 = temp;
    end
    
    z = (x1 + x2) / 2;
    fz = objective_function(z);
    dfz = first_derivative_function(z);
    funeval = funeval + 4;
    
    fprintf('#It\t\tx1\t\tx2\t\tf''(x1)\t\tf''(x2)\t\tf''(z)\t\tf(z)\n');
    
    while abs(dfz) >= epsilon
        fprintf('%d\t%.3f\t%.3f\t%10.3f\t%10.3f\t%10.3f\t%f\n', funeval, min(x1, x2), max(x1, x2), ...
            dfx1, dfx2, dfz, fz);
        
        if dfz < 0
            x1 = z;
        elseif dfz > 0
            x2 = z;
        end
        
        dfx1 = first_derivative_function(x1);
        dfx2 = first_derivative_function(x2);
        z = (x1 + x2) / 2;
        fz = objective_function(z);
        dfz = first_derivative_function(z);
        funeval = funeval + 4;
    end
    
    extrema = 'minima';
    if x1 > x2
        extrema = 'maxima';
    end
    
    fprintf('The %s lies between (%.3f, %.3f)\n', extrema, min(x1, x2), max(x1, x2));
    fprintf('Total number of function evaluations: %d\n', funeval);
    fprintf('============================================\n\n');
end

function val = objective_function(x)
    val = (2 * x - 5)^4 - (x^2 - 1)^3;
end

function val = first_derivative_function(x)
    val = 8 * (2 * x - 5)^3 - 6 * x * (x^2 - 1)^2;
end

