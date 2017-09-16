t = linspace(0,2*pi,1000);
x = sin(t);
y = zeros(length(t-1));
for i = 2:length(t)
    y(i) = (x(i) - x(i-1))/(t(i) - t(i-1));
end

hold on
plot(t,x,'b')
plot(t,y,'r')
hold off
