function f = Monoamine(t,y,v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Initialize
f    = zeros(9,1);

for i=1:length(v)
    
    tmp              = eval( v(i).fname );
    f(v(i).idx(:,1)) = f( v(i).idx(:,1) ) + v(i).idx(:,2)*tmp;
        
end
 
end

%--------------------------------------------------------------------------
% Auxillary functions
%--------------------------------------------------------------------------


function v = v_1Var(x,param)
% 1-variable reaction
v = param(2)*x/(param(1)+x);

end

function v = v_2Var(x,y,param)
% 2-variable reaction
v = param(3)*x*y/( (param(1)+x)*(param(2)+y) );

end

function v = v_4Var(x,y,u,v,param)
% 4-variable reaction

v = v_2Var(x,y,param(1:3)) - v_2Var(u,v,param(4:6));

end

function r = release(y)
% Release due to autoreceptors

if y < 0.000768
    r = 1.5 - (0.5/0.000768)*y;
elseif y < 0.0023
    r = 1 - (y-0.000768)*(0.6)/(0.0023-0.000768);
else
    r = 0.4;
end

end

function f = fire(t)
% Fire data

f = 1;

end


function v = v_TPH(x,y,z,param)
% Tryptophan hydroxylase

% Auto-receptor effect
auto  = 1.5 - (z^2) / ( (0.000768)^2 + z^2 );
denom = (param(1) + y + (y^2)/param(3) )*( param(2) + x );

% Reaction velocity
v = param(4)*x*y*auto./denom;

end


function v = v_TH(x,y,u,v,param)
% Tyrosine hydroxylase

% Auto-receptor effect
auto  = 0.5 + param(5) / ( 8*((v/0.002024)^4) + 1);
denom = ( x*y+param(1)*x+param(1)*param(2)*(1+u/param(3)) )*...
        ( 1 + y/param(4));

% Reaction velocity
v = param(5)*0.56*x*y*auto/denom;

end


function f = fluox(t)
% Fluox

f = 1;

end

