                        %% Project phase 2%%
   
clear all;
clc;

% Qu.1,n=5; Qu.2,n=3; Qu.3,n=4; Qu.4,n=6; Qu.5,n=2;

que=5; % Question number

% x0 = [-2;4]    % Intial guess vector
 %x0 = [-0;0;0]
%  x0 = [5;6;7;8]
  x0 = [-1;-3;1;-3;0] 
%  x0 = [-2;4;2;2;-4;5]


[Optimum_point, Function_value, Iteration] = project_phase_2(x0,que)


function [answer1, answer2, answer3]= project_phase_2(x0,que)
    % x0 = Initial Approximation %
    k = 1 ; % Iteration Counter
    e = 0.00001 ; % Accuracy
while true
    gradient_0 = grad(x0,que) ;  % Gradient at initial Point
if (magnitude(gradient_0) < e) % Termination condition 1 
    answer1 = x0 ;
    answer2= Function(x0,que);
    answer3 = k;
    break
end

direction = gradient_0 ; % Direction for unidirectional Search 
Alpha = project_1(x0,direction,que) ; % Value of Alpha by unidirectional Search using Project Phase 1
x1 = new_point(x0,Alpha,direction) ; % New point after substituting the value
if (((magnitude(x1 - x0))/(magnitude(x0))) < e) % Termination COndition 2
    answer1 = x1;
    answer2 = Function(x1,que);
    answer3 = k;
    break
end 
    k = k + 1;
    x0 = x1;
end
end

% ------- Multivariable Function ------------%
function fun_val = Function(x,que)
    fun_val=0;
    n=5;
if (que==1)
    %n=3;
for i=1:n
    fun_val=fun_val+(i*x(i))^2;
end
elseif (que==2)
    %n=3;
    for i=1:(n-1)
    fun_val= fun_val + 100*(x(i+1)-(x(i))^2)^2 + (x(i)-1)^2;
    end
elseif (que==3)
    %n=3;
    fun_val= (x(1)-1)^2;
    for i=2:n
    fun_val= fun_val + i*(2*(x(i))^2 - x(i-1))^2 ;
    end
elseif (que==4)
    %n=6;
    fun_val=(x(1)-1)^2;
    for i=2:n
      fun_val=fun_val + (x(i)-1)^2 - x(i)*x(i-1);
    end
elseif (que==5)
    %n=2;
    a=0;
    b=0;
    c=0;
    for i=1:n
        a= a + (x(i))^2;
        b= b + 0.5*i*x(i);
        c= c + 0.5*i*x(i);
    end
    fun_val= a + b^2 + c^4; 
else
    fun_val=0;
end 
end

%---------Function for gradient ------------%
function gradient = grad(x,que) 
    gradient = zeros(length(x),1);
    h = 0.001;
for i = 1:length(x)
    y = x;
    y(i) = y(i)+h;
    a = Function(y,que);
    y(i) = y(i)-2*h;
    b = Function(y,que);
    gradient(i) = (a - b)/(2*h);
end
end

%---Fuction creation for single variable optimization for gradient---%
function fun_val = Objective_Fun(y,x0,direction,que)
    l = length(x0);
for i = 1:l
    x0(i) = x0(i) - y*direction(i);
end
    fun_val = Function(x0,que);
end

%-------------Function for magnitude of a vector--------------%
function m = magnitude(gradient)
    magnitude_squeuare = gradient.*gradient;
    magnitude_squeuare_sum = sum(magnitude_squeuare);
    m = sqrt(magnitude_squeuare_sum);
end

%------Function to find new vector after finding the value of Alpha-------%
function new_points = new_point(x,Alpha,direction)
    l = length(x) ;
for i = 1:l
    x(i) = x(i) - Alpha*direction(i) ;
end
    new_points = x ;
end

%------------Function of single variable optimization----------%
%------------------Same as Project Phase 1 --------------------%

function Alpha = project_1(x,direction,que)
    % Exaustive search method %
    epsilon=0.001;
    a=-10;
    b=10;
    n=100;
    delta=(b-a)/n;
    x1=a;
    x2=x1+delta;
    x3=x2+delta;
    fd1= Objective_Fun(x1,x,direction,que);
    fd2= Objective_Fun(x2,x,direction,que);
    fd3= Objective_Fun(x3,x,direction,que);
while(x3<=b)
    if ((fd1 >= fd2 && fd2 <= fd3))
        %ALpha1=[x1 x3];
        break;
    else
        x1 = x2;
        x2 = x3;
        x3 = x2 + delta;
        fd1 = fd2;
        fd2 = fd3;
        fd3= Objective_Fun(x3,x,direction,que); 
    end

end
    %Alpha1=[x1 x3];

% Bisection method %

    x11=x1;
    x12=x3;
    h=0.001;
    % xd11 = x11 + h;
    % xd12 = x11 - h;
    % dfn1 =(Objective_Fun(xd11,x,direction,que)- Objective_Fun(xd12,x,direction,que))/(2*h);
    % xd21 = x12 + h;
    % xd22 = x12 - h;
    % dfn2 =(Objective_Fun(xd21,x,direction,que)- Objective_Fun(xd22,x,direction,que))/(2*h);

    z=(x11+x12)/2;
    zd21 = z + h;
    zd22 = z - h;
    dfz =(Objective_Fun(zd21,x,direction,que)- Objective_Fun(zd22,x,direction,que))/(2*h);

while(z >= a && z <= b)
    if (abs(dfz) < epsilon) 
        Alpha=z;
        break;
    elseif (dfz < 0)   
        x11 = z;
    elseif (dfz > 0) 
        x12 = z;
    end
    z=(x11+x12)/2;
    zd21 = z + h;
    zd22 = z - h;
    dfz =(Objective_Fun(zd21,x,direction,que)- Objective_Fun(zd22,x,direction,que))/(2*h);        
end
    Alpha=z;
end

