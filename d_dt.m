function [xd] = d_dt(x,dt)
% numerical differentiation with second order taylor series approximation. 
s=size(x);
if s(2)>s(1)
x = [x x(end)+(x(end)-x(end-1))+((x(end)-x(end-1))-(x(end-1)-x(end-2)))];
xd = diff(x)./dt;
elseif s(1)>s(2)
x = [x;x(end)+(x(end)-x(end-1))+((x(end)-x(end-1))-(x(end-1)-x(end-2)))];
xd = diff(x)./dt;
end
end

