%% HW 2 Problem 3
clear all; clc
%% Parameters from problem statement
a_0 = 26600;
i_0 = 1.10654;
e_0 = .74;
w_0 = 5*pi/180;
O_0 = 90*pi/180;
M_0 = 10*pi/180;

% J2 of Earth
J_2 = 0.00108;
% Equatorial radius of the Earth
R = 6370;
% Gravitational Parameter
mu = 3.986e05;

%time of propogation
t_end = 100*24*60*60;

%% Given above -> Initial Calculations
h_0 = sqrt(mu*a_0*(1 - e_0^2));

%% Convert Mean Anomaly to True Anamoly
%% Implementation of Laguerre-Conway: My own code from a previous class
n = 5;
E = M_0;
M = M_0;
e = e_0;

del = 1000000000;
count = 0;
while(del > .000001)
    F = E - e*sin(E) - M;
    F_1 = 1 - e*cos(E);
    F_2 = e*sin(E);
    E_old = E;
    %Laguerre-Conway
    E = E - (n * F)/ ...
        (F_1+sign(F_1)*abs((n-1)^2*(F_1)^2 - n*(n-1)*F*F_2)^.5);
    del = (E - E_old);
    count  = count + 1;
end

% True Anamoly from eccentric anamoly
temp = tan(E/2) * sqrt( (1 + e) / (1 - e) );
TA = 2 * atan(temp);
TA_0 = TA;

%% Perturbed Equations
a = a_0;
i = i_0;
w = w_0;
h = h_0;
u_0 = w_0 + TA_0;
O = O_0;

options = odeset(...
 'reltol', 1.e-10, ...
 'abstol', 1.e-10, ...
'initialstep', 10);

%% This method is similar to a method in Curtis Matlab Example 10_6.
%  Integrating the rates given the intial conditions of the problem
elements_0 = [h_0, e_0, TA_0, O_0, i_0, w_0];
y0 = elements_0';
t_array = 0:10:t_end;
[t, y] = ode89(@rates, t_array, y0, options);

h  = y(:, 1); 
e  = y(:, 2);
TA = y(:, 3);
O  = y(:, 4);
i  = y(:, 5); 
w  = y(:, 6);

a = h.^2./mu .* 1 ./ (1 - e.^2); 

%% Plotting Orbital Elements
figure
subplot(6, 1, 1)
hold on
plot(t_array/(3600*24), a); title('Semimajor Axis (a)')
ylabel('Semimajor Axis [km]'); xlabel(['Time [days]']); grid minor;
hold off

subplot(6, 1, 2)
hold on
plot(t_array/(3600*24), i*180/pi); title('Inclination (i)')
ylabel('Inclination [deg]'); xlabel(['Time [days]']); grid minor;
hold off

subplot(6, 1, 3)
hold on
plot(t_array/(3600*24), e); title('Eccentricity (e)')
ylabel('Eccentricity'); xlabel(['Time [days]']); grid minor;
hold off

subplot(6, 1, 4)
hold on
plot(t_array/(3600*24), w*180/pi); title('Argument of Periapsis (w)')
ylabel('Argument of Periapsis [deg]'); xlabel(['Time [days]']); grid minor;
hold off

subplot(6, 1, 5)
hold on
plot(t_array/(3600*24), O*180/pi); title('Right Acension (O)')
ylabel('Right Acension [deg]'); xlabel(['Time [days]']); grid minor;
hold off

subplot(6, 1, 6)
hold on
plot(t_array/(3600*24), h); title('Angular Momentum (h)')
ylabel('Angular Momentum [km^2/sec]'); xlabel(['Time [days]']); grid minor;
hold off

%% Subfunction that generates rates
function dfdt = rates(t, f)
    % Earth Constants
    mu = 3.986e05;
    J_2 = 0.00108;
    R = 6370;

    h  = f(1);
    e  = f(2);
    TA = f(3);
    O  = f(4);
    i  = f(5);
    w  = f(6);

    % Curtis pg 511, Get r, u from parameters
    r = h^2/(mu * (1 + e*cos(TA)));
    u = w + TA;

    % Gravitational Perturbation in the RSW frame, Curtis 10.88
    p_r = (-3/2) * (J_2*mu*R^2) / (r^4) * (1 - 3*sin(i)^2 * sin(u)^2); 
    p_s = (-3/2) * (J_2*mu*R^2) / (r^4) * sin(i)^2 * sin(2*u);
    p_w = (-3/2) * (J_2*mu*R^2) / (r^4) * sin(2*i) * sin(u);
    
    % Gauss Planetary Equations, Curtis 10.84
    hdot = r*p_s;
    edot = h/mu * sin(TA) * p_r + 1/(mu*h) * ((h^2 + mu*r) * cos(TA) + mu*e*r) * p_s;
    TAdot = h/(r^2) + 1/(e*h) * (h^2/mu * cos(TA) * p_r - (r + h^2/mu) * sin(TA) * p_s);
    Odot = r/(h*sin(i)) * sin(u) *p_w;
    idot = r/h * cos(u) * p_w;
    wdot = -1/(e*h) * (h^2/mu * cos(TA) *p_r - (r + h^2/mu) * sin(TA) * p_s) - ...
        ( r*sin(u) ) / ( h * tan(i) ) * p_w;

   
    dfdt = [hdot, edot, TAdot, Odot, idot, wdot]';
end
