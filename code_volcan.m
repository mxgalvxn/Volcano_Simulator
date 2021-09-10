%simulacion volcan
clc;
clear;
set(gcf, 'Position', get(0,'Screensize')); warning('off','all'); %Pantalla Completa
hold on; 
grid on;
plot(-1,-1,'-y',-2,-2,'MarkerSize',10); legend('Tiro Parabolico');
set(0,'DefaultLegendAutoUpdate','off');
axis([-8000 8000 0 5000]);
delete(findall(gcf,'type','text'));

title(('Simulacion de proyectil en tiro parabolico'));
xlabel('Distancia [m]');
ylabel('Altura [m]');

x0 = 0; %Altura inicial en x (m)
y0 = 2000; %Altura inicial en y (m)
g = 9.81;

a=-750:.1:750; %rango volcan
b=2000*cos(.0021*a);%función volcan
plot(a,b,'r')%dibujo volcan
hold on
posx = 3100;
posy = 3000;
for i= 1:1:3
v0= 10 + (200)*rand(1);
ang =  20 + (160)*rand(1);
angrad = deg2rad(ang);

deltat = 0.1;
cdrag = .011;
dens_aire = 1.2;
dens_roca = 3120;
masa = 2000;
vol = masa / dens_roca;
radio =((3*vol)/4*pi)^(1/3);
area = (pi * radio^2) / 4;
b = 0.5 * dens_aire * cdrag * area;
t = 0:deltat:(deltat)^2;
vx=cos(angrad)*(v0);
vy=sin(angrad)*(v0);

%x
x_=vx*t;
y_ = vy*t+0.5*(-g)*(t.^2)+y0;

dvx = -abs(vx) / vx * b / masa * vy^2;
x_menos = x_(1) - vx * deltat - (dvx * deltat^2);
x = zeros(1,100);
x(1) = x_menos;
x(2) = x_(1);
%y
dvy = g - abs(vy) / vy * b / masa * vy^2;
y_menos = y_(1) - vy * deltat - dvx * deltat^2;
y = zeros(1,100);
y(1) = y_menos;
y(2) = y_(1);

for z = 3:1:500
x(z) = (2 .* x(z-1) - x(z-2)) - ((x(z-1) - x(z-2)) / (deltat))^2 * (b / masa) * (deltat^2)
y(z) = (2 .* y(z-1) - y(z-2)) - ((y(z-1) - y(z-2)) / (deltat))^2 * (b / masa) * (deltat^2) - (g * deltat^2)
if y(z)<=0
    break
end

end
poly = [0.5*-g, vy, y0];
sol = roots(poly);
if(sol(1)>0)
    tvuelo = sol(1);
else
    tvuelo = sol(2);
end

comet(x,y)

hold on 
text(3000,posy-150,['Altura maxima = ' num2str(max(y))]);
text(3000,posy-250,['Distancia maxima = ' num2str(max(x))]);
text(3000,posy-350,['Angulo = ' num2str(ang) '°'] );
text(3000,posy-450,['Velocidad inicial = ' num2str(v0)]);
text(3000,posy-550,['Tiempo = ' num2str(tvuelo)]);
posy = posy - 600;
t4 = text(x(length(x)),y(length(y)), ['Roca ' num2str(i)]);
end



text(3000,3200,['Altura inicial =' num2str(y0)]);
