clear;
x1 = [-30 -20 -10.0 0.0 10.0 20.0 35.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 110.0 ];
y1 = [0.0 -10.0 -25.0 -60.0 -90.0 -105.0 -80.0 -55.0 -30.0 -15.0 0.0 10.0 30.0 40.0 65.0];
a = 176.0914
b = -1722.096
c = -197.758
d = 1956.205
for I = 1:15
   y2(I) = a + c / (1 + exp(-b*(x1(I) - d)))
end
plot(x1,y1,x1,y2);
