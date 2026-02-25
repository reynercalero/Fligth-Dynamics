clear all; clc;

% Cambiar indices segun se requiera
caso_vuelo = 'A';

switch caso_vuelo
    case 'A' % vuelo recto y nivelado
        phi = 0;
        theta = 5;
        psi = 0;
        Vb = [650; 0; 100];
    case 'B' % Climb
        phi = 0;
        theta = 10;
        psi = 0;
        Vb = [98; 0; 5];
    case 'C' % Aircraft Turn
        phi = 30;
        theta = 5;
        psi = 45;
        Vb = [100; 2; 0];
end

R = matriz_rotacion(phi, theta, psi);
Vned = R * Vb;
[alpha, beta, gamma] =angulos_aeros(Vb, theta);

fprintf('=== Resultados caso de vuelo %s ===\n', caso_vuelo);
fprintf('Angulos Euler [deg]: Roll: %.1f, Pitch: %.1f, Yaw: %.1f\n', phi, theta, psi);
fprintf('Velocidad body [u,v,w]: [%.1f, %.1f, %.1f]\n', Vb(1), Vb(2), Vb(3));
fprintf('------------------------------------\n');
fprintf('Velocidad NED [N,E,D]:  [%.2f, %.2f, %.2f] Km/h \n', Vned(1), Vned(2), Vned(3));
fprintf('Angulos aero: Alpha: %.2f°, Beta: %.2f°, Gamma: %.2f°\n', alpha, beta, gamma);
visualizador(R, Vb, caso_vuelo)

