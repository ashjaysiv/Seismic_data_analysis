function [A, T] = responsespectrum(accel, ee,dt)

% ee - damping in % - 5 is recommended
% y - gamma in newmark's method - 0.5 is recommended
% b - beta in newmark's method - 0.25 is recommended
% td - time till which you want graph to be plotted

Tn = 10; % time period till which u want response spectrum
y = 0.5;
b = 0.25;
uo = 0;
vo = 0;
m = 1;
z = ee/100;
% dt = 0.005;
% dt = 0.02;

na = length(accel);
nl = 2*na;

T = [0.01;0.015;0.02;0.03;0.04;0.05;0.06;0.075;0.09;0.1;0.15;0.2;0.3;0.4;0.5;0.6;0.7;0.75;0.8;0.9;1;1.2;1.5;2;2.5;3;4;5;6;7.5;8;9;10];
accel = [accel;zeros(nl-na,1)];
p = -m*accel;

A = zeros(length(T),1);% acclelration response spectrum - total accelaration
V = zeros(length(T),1);% velocity response spectrum - relative velocity
D = zeros(length(T),1);% displacement response spectrum - relative displacement
for j = 1:length(T)
    
    fn = 1/T(j);
    wn = 2*pi*fn;
    k = m*wn^2;
    c = 2*m*wn*z;
    
    u = zeros(nl,1);
    v = zeros(nl,1);
    ac = zeros(nl,1);
    
    u(1) = uo;
    v(1) = vo;
    ac(1) = (p(1)-c*vo-k*uo)/m;
    
    kf = k + y*c/(b*dt) + m/(b*dt^2);
    a = m/(b*dt) + y*c/b;
    b2 = m/(2*b) + dt*(y/(2*b) - 1)*c;
    
    for i = 1:nl-1
        p1=p(i);
        p2=p(i+1);
        dpf = (p2 - p1) + a*v(i) + b2*ac(i);
        du = dpf/kf;
        dv = y/(b*dt)*du - (y/b)*v(i) + dt*(1 - y/(2*b))*ac(i);
        da = du/(b*dt^2) - v(i)/(b*dt) - ac(i)/(2*b);
        u(i+1) = u(i) + du;
        v(i+1) = v(i) + dv;
        ac(i+1) = ac(i) + da;
    end
    
    asd = ac + accel;
    A(j) = max(abs(asd));
    V(j) = max(abs(v));
    D(j) = max(abs(u));
end

A = [max(abs(accel(:)));A];
V = [0;V];
D = [0;D];

PSV = (2*pi./T).*D(2:end); % pseudo spectral velocity
PSV = [PSV(1); PSV];
PSA = ((2*pi./T).^2).*D(2:end); % pseudo spectral accleration
PSA = [PSA(1); PSA];

T = [0;T];

% figure
% plot(T, A)
% xlabel 'Time Period (seconds)'
% ylabel 'Spectral Accelaration'
% 
% figure
% plot(T, V)
% xlabel 'Time Period (seconds)'
% ylabel 'Spectral Velocity'
% 
% figure
% plot(T, D)
% xlabel 'Time Period (seconds)'
% ylabel 'Spectral Displacement'
% 
% figure
% plot(T, PSV)
% xlabel 'Time Period (seconds)'
% ylabel 'Pseudo Spectral Velocity'
% 
% figure
% plot(T, PSA)
% xlabel 'Time Period (seconds)'
% ylabel 'Pseudo Spectral Accelaration'

end
