%%HW 2 Problem 2
clear
J_2 = 0.00196; R = 3390; mu = 4.282e04;

%Problem 1
%Want to orbit the Earth three times a day
T_p = 24*60*60 + 39*60 + 35; 
n   = 2*pi/T_p;
a   = ((T_p/(2*pi))^2*mu)^(1/3);
i   = acos(sqrt(.2));
e_  = [0:.001:1];
e   = 1 - (R+400)/a;

%O_dot_ as a function of e_
O_dot_ = -3/2 .* n .* J_2 .* (R/a)^2 .* cos(acos(1/sqrt(5))) ./ ((1 - e_.^2) .^ 2);

figure
hold on
plot(e_, O_dot_); xlabel('Eccentricity'); ylabel('$\dot{\overline{\Omega}}$ (rad/sec)', ...
    'Interpreter','Latex');
title('$\dot{\overline{\Omega}}$(e)', 'Interpreter','Latex'); grid minor; ylim([-.0001 .0001])
hold off

%Solution of O_dot1 given the choice of a, e, i
O_dot1 = -3/2 .* n .* J_2 .* (R/a)^2 .* 1/sqrt(5) ./ ((1 - e.^2) .^ 2);