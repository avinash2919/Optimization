                        %% Project phase 3 %%
   
clear;
clc;

  que = 4; % Question number

  x0 = [0;0]    % Intial guess vector

% x0 = [-2;4;2;2;-4;5;2;5]


[Optimum_point, Function_value, Iteration] = project_phase_3(x0,que)

function [answer1, answer2, answer3]= project_phase_3(x0,que)
    % x0 = Initial Approximation %
    t = 0 ; % Iteration Counter
    k = 1 ; % Iteration Counter
    c = 10;
    r = 0.1;
    e = 0.0001 ; % Accuracy

while true
    gradient_0 = grad(x0,que,r);   % Gradient at initial Point
if (magnitude(gradient_0) < e) % Termination condition 1 
        x1=x0;
%     answer1 = x0 ;
%     answer2 = Function(x0,que,r);
%     answer3 = k;
    break
end

direction = gradient_0;  % Direction for unidirectional Search 
Alpha = project_1(x0,direction,que,r);  % Value of Alpha by unidirectional Search using Project Phase 1
x1 = new_point(x0,Alpha,direction);  % New point after substituting the value
if (((magnitude(x1 - x0))/(magnitude(x0))) < e) % Termination COndition 2
%     answer1 = x1;
%     answer2 = Function(x1,que,r);
%     answer3 = k;
    break
end 
    k = k + 1;
    x0 = x1;
end
y = Function(x1,que,r);

while true
    r = c*r;
    while true
        gradient_0 = grad(x0,que,r) ;  % Gradient at initial Point
        if (magnitude(gradient_0) < e) % Termination condition 1
              x1=x0;
%             answer1 = x0 ;
%             answer2 = Function(x0,que,r);
%             answer3 = k;
        break
        end
    
        direction = gradient_0 ; % Direction for unidirectional Search 
        Alpha = project_1(x0,direction,que,r) ; % Value of Alpha by unidirectional Search using Project Phase 1
        x1 = new_point(x0,Alpha,direction) ; % New point after substituting the value
        if (((magnitude(x1 - x0))/(magnitude(x0))) < e) % Termination COndition 2
                x1;
%             answer1 = x1;
%             answer2 = Function(x1,que,r);
%             answer3 = k;
        break
        end 
        k = k + 1;
        x0 = x1;
    end
        z=Function(x1,que,r);
        if (abs(z-y)<e)
            answer1 = x1;
            answer2 = z; %Function(x1,que,r);
            answer3 = t+k;
        break
        end 
        t = t + 1;
        y=z;
        x0=x1;
        answer1 = x1;
        answer2 = z; %Function(x1,que,r);
        answer3 = t+k;
end
end

%---------Function for bracket ------------%
function bracket = brak(g) %constraint of >= type.
    if(g<0)
    bracket=g;
    else
    bracket=0;    
    end    
end
% ------- Multivariable Function ------------%
function fun_val = Function(x,que,r)

if (que==1)

    g1=((x(1)-5)^2+(x(2)-5)^2-100);
%      if(g1<0)
%      f1=g1;
%      else
%      f1=0;    
%      end 
     g2=(-(x(1)-6)^2-(x(2)-5)^2+82.81);
%      if(g2<0)
%      f2=g2;
%      else
%      f2=0;    
%      end 
    
    j=r*((brak(g1))^2+r*((brak(g2))^2));
    fun_val=((x(1)-10)^3+(x(2)-20)^3)+j;
%   fun_val=((x(1)-10)^3+(x(2)-20)^3)+r*((f1)^2+(f2)^2);
   

elseif (que==4)
    g1=((x(1)-5)^2+(x(2))^2-26);
    fun_val=((x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2)+r*((brak(g1))^2);
elseif (que==2)
    g1=(-x(1)^2+x(2)-1);
    g2=-1+x(1)-(x(2)-4)^2;
    fun_val=((sin(2*pi*x(1)))^3)*(sin(2*pi*x(2)))/((x(1)^3)*(x(1)+x(2)))+r*((brak(g1))^2+r*((brak(g2))^2));
elseif (que==3)
    g1=double((x(1)-5)^2+(x(2)-5)^2-100);
    g2=double(-(x(1)-6)^2-(x(2)-5)^2+82.81);
    g3=double((x(1)-5)^2+(x(2)-5)^2-100);
    g4=double(-(x(1)-6)^2-(x(2)-5)^2+82.81);
    g5=double((x(1)-5)^2+(x(2)-5)^2-100);
    g6=double(-(x(1)-6)^2-(x(2)-5)^2+82.81);
    fun_val=((x(1)-10)^3+(x(2)-20)^3)+r*((brak(g1))^2+r*((brak(g2))^2)+(brak(g3))^2+r*((brak(g4))^2)+(brak(g5))^2+r*((brak(g6))^2));
else
    fun_val=0;
end 
end

%---------Function for gradient ------------%
function gradient = grad(x,que,r) 
    gradient = zeros(length(x),1);
    h = 0.001;
for i = 1:length(x)
    y = x;
    y(i) = y(i)+h;
    a = Function(y,que,r);
    y(i) = y(i)-2*h;
    b = Function(y,que,r);
    gradient(i) = (a - b)/(2*h);
end
end

%---Fuction creation for single variable optimization for gradient---%
function fun_val = Objective_Fun(y,x0,direction,que,r)
    l = length(x0);
for i = 1:l
    x0(i) = x0(i) - y*direction(i);
end
    fun_val = Function(x0,que,r);
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

function Alpha = project_1(x,direction,que,r)
    % Exaustive search method %
    epsilon=0.001;
    a=-10;
    b=10;
    n=100;
    delta=(b-a)/n;
    x1=a;
    x2=x1+delta;
    x3=x2+delta;
    fd1= Objective_Fun(x1,x,direction,que,r);
    fd2= Objective_Fun(x2,x,direction,que,r);
    fd3= Objective_Fun(x3,x,direction,que,r);
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
        fd3= Objective_Fun(x3,x,direction,que,r); 
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
    dfz =(Objective_Fun(zd21,x,direction,que,r)- Objective_Fun(zd22,x,direction,que,r))/(2*h);

while(z >= a && z <= b)
    if (abs(dfz) < epsilon) 
        %Alpha=z;
        break;
    elseif (dfz < 0)   
        x11 = z;
    elseif (dfz > 0) 
        x12 = z;
    end
    z=(x11+x12)/2;
    zd21 = z + h;
    zd22 = z - h;
    dfz =(Objective_Fun(zd21,x,direction,que,r)- Objective_Fun(zd22,x,direction,que,r))/(2*h);        
end
    Alpha=z;
end

