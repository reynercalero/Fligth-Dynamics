[alpha, beta, gamma, t] = angulos_aeros();
[t,phi,theta,psi,phi_dt,theta_dt,psi_dt]=attitude();
[t,px,py,pz,Vx,Vy,Vz,Xac,Yac,Zac,a_ned,a_real]=posicion();

##figure(1);
##plot(t, alpha, 'r', 'DisplayName', 'Alpha'); hold on;
##plot(t, beta, 'b', 'DisplayName', 'Beta');
##plot(t, gamma, 'g', 'DisplayName', 'Gamma');
##grid on;
##xlabel('Tiempo (s)');
##ylabel('Ángulos (deg)');
##legend()

figure(2);
subplot(3,1,1);
plot(t,rad2deg(phi),'r','lineWidth',1.2);
grid on;
ylabel('Roll \phi (deg)');
title('Actitud de la Aeronave (Ángulos de Euler)');
subplot(3,1,2);
plot(t, rad2deg(theta), 'g', 'LineWidth', 1.2);
grid on;
ylabel('Pitch \theta (deg)');
subplot(3,1,3);
plot(t, rad2deg(psi), 'b', 'LineWidth', 1.2);
grid on;
ylabel('Yaw \psi (deg)');
xlabel('Tiempo (s)');

figure(3);
plot3(py, px, pz, 'b', 'LineWidth', 2);
grid on;
axis equal;
xlabel('East [m]');
ylabel('North [m]');
zlabel('Altura (Altitude) [m]');
title('Recorrido Espacial en marco NED');

hold on;
plot3(py(1), px(1), pz(1), 'go', 'MarkerFaceColor', 'g');
plot3(py(end), px(end), pz(end), 'ro', 'MarkerFaceColor', 'r');
legend('Trayectoria', 'Inicio', 'Fin');

figure(4);
subplot(3,1,1);
plot(t, px, 'r','lineWidth',1.2);
grid on;
ylabel('North [m]');
title('Evolución Temporal de la Posición');

subplot(3,1,2);
plot(t, py, 'g','lineWidth',1.2);
grid on;
ylabel('East [m]');

subplot(3,1,3);
plot(t, pz, 'b','lineWidth',1.2);
grid on;
ylabel('Altura [m]');
xlabel('Tiempo [s]');
