function RCAM()
    clear; clc; close all;
    fprintf('  RCAM — Research Civil Aircraft Model\n');
    pkg load optim;

    % --- Inicialización Caso 1, 2 y 3 ---
    X0 = [85.0; 0.0; 0.0; ...   % u, v, w   [m/s]
           0.0; 0.0; 0.0; ...   % p, q, r   [rad/s]
           0.0; 0.1; 0.0];      % phi, theta, psi [rad]
    U0 = [0.0; -0.1; 0.0; 0.08; 0.08];

    % [CASO 1]
    fprintf('\n[CASO 1] \n');
    simulate(X0, U0, 180.0, 0.1, [], [], 'Caso 1');

    % [CASO 2] Deflexión de alerón
    fprintf('\n[CASO 2] Deflexión alerón +5 deg (t = 30..32 s)\n');
    simulate(X0, U0, 180.0, 0.1, 5.0, [], 'Caso 2 - Deflexion aleron +5 deg');

    % [CASO 4] Fallo motor 2
    fprintf('\n[CASO 3] Fallo motor 2\n');
    simulate(X0, U0, 180.0, 0.1, [], 2, 'Caso 3 - Fallo motor 2');

    % [CASO 5] PSO — Trim Va=78 m/s, psi=45 deg
    Va_obj = 78.0;
    psi_obj = pi / 4; % 45 deg

    fprintf('\n[CASO 4] PSO — Trim Va=78 m/s, psi=45 deg \n');
    [U_pso, J_pso, X_pso, hist_J] = pso_trim(Va_obj, psi_obj, 50, 500);

    % Punto inicial (X0 para optimización local: dim=7)
    alpha_pso = atan2(X_pso(3), X_pso(1));
    x0_refine = [U_pso(1); U_pso(2); U_pso(3); U_pso(4); U_pso(5); alpha_pso; X_pso(8)];

    options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxIter', 50000, 'MaxFunEvals', 100000);

    obj_refine_func = @(part) obj_refine_local(part, Va_obj, psi_obj);

    [x_ref, J_refine] = fminsearch(obj_refine_func, x0_refine, options);

    if J_refine < J_pso
        fprintf('  Refinamiento MEJORÓ: J_PSO = %.6f  ->  J_refine = %.6f\n', J_pso, J_refine);
        [X_trim, U_trim] = part_to_state_ctrl(x_ref, Va_obj, psi_obj);
        J_trim = J_refine;
    else
        fprintf('  Refinamiento no mejoró (J_refine = %.6f >= J_PSO = %.6f)\n', J_refine, J_pso);
        X_trim = X_pso; U_trim = U_pso; J_trim = J_pso;
    end

    imprimir_resultados_pso(U_trim, J_trim, X_trim);

    figure('Name', 'Convergencia PSO');
    semilogy(hist_J, 'b', 'LineWidth', 1.4);
    xlabel('Iteración'); ylabel('J (costo)');
    title('Convergencia PSO — Trim RCAM');
    grid on;

end


% funciones


function Xd = xdot(X, U)
    m = 120000.0; g = 9.81; S = 260.0; mac = 6.6; b = 44.8;
    S_t = 64.0; l_t = 24.8; rho = 1.225;
    alpha_0 = -11.5 * pi/180; n_lift = 5.5; deps_dalpha = 0.25;

    r_apt1 = [0.0;  7.94; 1.9]; r_apt2 = [0.0; -7.94; 1.9];
    r_cg = [0.23*mac; 0.0; 0.1*mac]; r_ac = [0.12*mac; 0.0; 0.0];

    Ib = m * [40.07,  0.0,   2.098; ...
               0.0,  64.0,   0.0; ...
               2.098, 0.0,  99.92];

    % Desempaquetado de estados y controles
    V_b = X(1:3); u_b = X(1); v_b = X(2); w_b = X(3);
    omega = X(4:6); p = X(4); q = X(5); r = X(6);
    phi = X(7); theta = X(8);

    da = max(min(U(1),  25*pi/180), -25*pi/180);
    de = max(min(U(2),  10*pi/180), -25*pi/180);
    dr = max(min(U(3),  30*pi/180), -30*pi/180);
    th1 = U(4); th2 = U(5);

    % Variables Intermedias
    Va = max(norm(V_b), 1e-6);
    alpha = atan2(w_b, u_b);
    beta = asin(max(min(v_b/Va, 1.0), -1.0));
    Q = 0.5 * rho * Va^2;

    % Coeficientes de fuerza aerodinámica
    if alpha <= 14.5 * pi/180
        CL_wb = n_lift * (alpha - alpha_0);
    else
        CL_wb = 15.212 - 155.2*alpha + 609.2*alpha^2 - 768.5*alpha^3;
    end

    epsilon = deps_dalpha * (alpha - alpha_0);
    alpha_t = alpha - epsilon + de + 1.3*q*l_t/Va;
    CL_t = (S_t/S) * 3.1 * alpha_t;

    CL = CL_wb + CL_t;
    CD = 0.13 + 0.07*(n_lift*alpha + 0.654)^2;
    CY = -1.6*beta + 0.24*dr;

    % Rotación Estabilidad -> Cuerpo
    F_As = Q * S * [-CD; CY; -CL];
    ca = cos(alpha); sa = sin(alpha);
    C_bs = [ ca, 0.0, -sa; ...
            0.0, 1.0, 0.0; ...
             sa, 0.0,  ca];
    F_Ab = C_bs * F_As;

    % Momentos Aerodinámicos en AC
    n_vec = [-1.4 * beta; ...
             -0.59 - (3.1*S_t*l_t)/(S*mac) * (alpha - epsilon); ...
             (1.0 - alpha*180.0/(pi*15.0)) * beta];

    dCm_domega = (mac/Va) * [-11.0,  0.0,   5.0; ...
                               0.0, -4.03*S_t*l_t^2/(S*mac^2), 0.0; ...
                               1.7,  0.0, -11.5];

    dCm_du = [-0.6,  0.0,  0.22; ...
               0.0, -3.1*(S_t*l_t)/(S*mac), 0.0; ...
               0.0,  0.0, -0.63];

    C_Mac = n_vec + dCm_domega * omega + dCm_du * [da; de; dr];
    M_Aac = C_Mac * Q * S * mac;

    % Momento respecto al CG
    M_Acg = M_Aac + cross(F_Ab, r_cg - r_ac);

    % Propulsión
    F1 = min(th1 * m * g, 0.175 * m * g);
    F2 = min(th2 * m * g, 0.175* m * g);
    F_Eb = [F1 + F2; 0.0; 0.0];
    M_Ecg = cross(r_apt1 - r_cg, [F1; 0.0; 0.0]) + cross(r_apt2 - r_cg, [F2; 0.0; 0.0]);

    % Gravedad
    F_gb = m * [-g * sin(theta); ...
                 g * cos(theta) * sin(phi); ...
                 g * cos(theta) * cos(phi)];

    % Ecuaciones diferenciales dinámicas (10 pasos)
    F_total = F_Ab + F_Eb + F_gb;
    M_total = M_Acg + M_Ecg;

    Xd_123 = F_total/m - cross(omega, V_b);
    Xd_456 = Ib \ (M_total - cross(omega, Ib * omega));

    % Cinemática de Euler
    sphi = sin(phi); cphi = cos(phi);
    stheta = sin(theta); ctheta = cos(theta); ttheta = tan(theta);

    H = [1.0,  sphi*ttheta,  cphi*ttheta; ...
         0.0,  cphi,        -sphi; ...
         0.0,  sphi/ctheta,  cphi/ctheta];
    Xd_789 = H * omega;

    Xd = [Xd_123; Xd_456; Xd_789];
end

function [t, X_hist] = simulate(X0, U0, t_total, dt, aileron_deg, engine_off, title_str)
    N = round(t_total / dt) + 1;
    t = linspace(0.0, t_total, N)';
    X_hist = zeros(N, 9);
    X_hist(1, :) = X0';
    X = X0;

    fprintf('  Iniciando: %s  (%d pasos, dt=%0.2fs)\n', title_str, N-1, dt);

    for k = 2:N
        t_k = t(k-1);
        U = U0;

        % defleccion aleron
        if ~isempty(aileron_deg) && (t_k >= 30.0 && t_k <= 32.0)
            U(1) = U(1) + aileron_deg * pi/180.0;
        end

        % Fallo de motor
        if ~isempty(engine_off)
            if engine_off == 1, U(4) = 0.0; end
            if engine_off == 2, U(5) = 0.0; end
        end

        % integracion
        Xd = xdot(X, U);
        X = X + dt * Xd;
        X_hist(k, :) = X';
    end
    fprintf('  Finalizado.\n');
    plot_results(t, X_hist, title_str);
end

function plot_results(t, X, title_str)
    figure('Name', title_str);

    subplot(3,3,1); plot(t, X(:,1), 'b'); ylabel('u (m/s)'); grid on; title([title_str '']);
    subplot(3,3,4); plot(t, X(:,2), 'b'); ylabel('v (m/s)'); grid on;
    subplot(3,3,7); plot(t, X(:,3), 'b'); ylabel('w (m/s)'); xlabel('Time (s)'); grid on;

    subplot(3,3,2); plot(t, X(:,4), 'b'); ylabel('p (rad/s)'); grid on; title('');
    subplot(3,3,5); plot(t, X(:,5), 'b'); ylabel('q (rad/s)'); grid on;
    subplot(3,3,8); plot(t, X(:,6), 'b'); ylabel('r (rad/s)'); xlabel('Time (s)'); grid on;

    subplot(3,3,3); plot(t, X(:,7), 'b'); ylabel('phi (rad)'); grid on; title('');
    subplot(3,3,6); plot(t, X(:,8), 'b'); ylabel('theta (rad)'); grid on;
    subplot(3,3,9); plot(t, X(:,9), 'b'); ylabel('psi (rad)'); xlabel('Time (s)'); grid on;
end

# FUNCION DE COSTO

function J = cost_function(Va_ideal, psi_ideal, X_eval, U)
    Xd = xdot(X_eval, U);
    Va_actual = norm(X_eval(1:3));
    alpha_actual = atan2(X_eval(3), X_eval(1));
    gamma_actual = X_eval(8) - alpha_actual;

    % Normalización y residuos adimensionales
    e_lineal = Xd(1:3) / Va_ideal;
    e_angular = Xd(4:6) / 0.1;
    e_euler = Xd(7:9) / 0.1;

    e_Va = (Va_ideal - Va_actual) / Va_ideal;
    e_gamma = gamma_actual / (5*pi/180);
    e_phi = X_eval(7) / (5*pi/180);
    e_psi = (psi_ideal - X_eval(9)) / (pi/4);

    errores = [e_lineal; e_angular; e_euler; e_Va; e_gamma; e_phi; e_psi];
    pesos = [200; 200; 200; ... % Xdot lineal
             100; 150; 100; ... % Xdot angular
              50;  50;  50; ... % Xdot Euler
             300;           ... % e_Va (Peso dominante)
             100;           ... % e_gamma
              80;           ... % e_phi
              80];              % e_psi

    J = sum(pesos .* (errores.^2));
end

function J = obj_refine_local(part, Va_ideal, psi_ideal)
    % Envuelve la cost_function para fminsearch
    [X, U] = part_to_state_ctrl(part, Va_ideal, psi_ideal);
    J = cost_function(Va_ideal, psi_ideal, X, U);
end

function [X, U] = part_to_state_ctrl(part, Va_ideal, psi_ideal)
    % Transforma el vector continuo de optimización de dim 7 a estado(9) y control(5)
    U = [part(1); part(2); part(3); part(4); part(5)];
    alpha = part(6); theta = part(7);
    u_b = Va_ideal * cos(alpha);
    w_b = Va_ideal * sin(alpha);
    X = [u_b; 0.0; w_b; 0.0; 0.0; 0.0; 0.0; theta; psi_ideal];
end
% OPTIMIZADOR PSO ADAPTATIVO
function [U_opt, gbest_cost, X_trim, historial_J] = pso_trim(Va_ideal, psi_ideal, n_particles, n_iter)
    d2r = pi / 180;

    % Límites físicos inferiores y superiores de búsqueda
    lb = [-25*d2r, -25*d2r, -30*d2r, 0.01, 0.01,  0*d2r, -5*d2r];
    ub = [ 25*d2r,  10*d2r,  30*d2r, 0.175, 0.175, 15*d2r, 15*d2r];
    d = 7;

    % Semilla fija para reproducibilidad
    rand('seed', 42);

    % Inicialización de partículas
    pos = repmat(lb, n_particles, 1) + rand(n_particles, d) .* repmat((ub - lb), n_particles, 1);
    v_max = 0.20 * (ub - lb);
    vel = repmat(-v_max, n_particles, 1) + rand(n_particles, d) .* repmat(2*v_max, n_particles, 1);

    % Evaluación Inicial
    cost = zeros(n_particles, 1);
    for i = 1:n_particles
        [X, U] = part_to_state_ctrl(pos(i,:), Va_ideal, psi_ideal);
        cost(i) = cost_function(Va_ideal, psi_ideal, X, U);
    end

    pbest_pos = pos;
    pbest_cost = cost;

    [gbest_cost, gbest_idx] = min(pbest_cost);
    gbest_pos = pbest_pos(gbest_idx, :);

    historial_J = zeros(n_iter, 1);
    fprintf('  %6s  %12s  %7s  %7s  %7s\n', 'Iter', 'J_best', 'w', 'c1', 'c2');

    % Bucle principal del PSO Adaptativo
    for t = 1:n_iter
        % Dinámica adaptativa de los hiperparámetros
        w = 0.4 * (t - n_iter) / (n_iter^2) + 0.4;
        c1 = -3.0 * (t / n_iter) + 3.5;
        c2 = 3.0 * (t / n_iter) + 0.5;

        r1 = rand(n_particles, d);
        r2 = rand(n_particles, d);

        gbest_pos_mat = repmat(gbest_pos, n_particles, 1);
        v_max_mat = repmat(v_max, n_particles, 1);
        ub_mat = repmat(ub, n_particles, 1);
        lb_mat = repmat(lb, n_particles, 1);

        % Actualización de velocidades y posiciones
        vel = w * vel + c1 * r1 .* (pbest_pos - pos) + c2 * r2 .* (gbest_pos_mat - pos);
        vel = max(min(vel, v_max_mat), -v_max_mat);

        pos = pos + vel;
        pos = max(min(pos, ub_mat), lb_mat);

        % Evaluación
        for i = 1:n_particles
            [X, U] = part_to_state_ctrl(pos(i,:), Va_ideal, psi_ideal);
            cost(i) = cost_function(Va_ideal, psi_ideal, X, U);
        end

        % Actualización mejores personales
        mejora = cost < pbest_cost;
        pbest_pos(mejora, :) = pos(mejora, :);
        pbest_cost(mejora) = cost(mejora);

        % Actualización mejor global
        [current_best_cost, idx] = min(pbest_cost);
        if current_best_cost < gbest_cost
            gbest_cost = current_best_cost;
            gbest_pos = pbest_pos(idx, :);
        end

        historial_J(t) = gbest_cost;

        if mod(t, 50) == 0 || t == 1
            fprintf('  %6d  %12.4f  %0.4f  %0.4f  %0.4f\n', t, gbest_cost, w, c1, c2);
        end
    end

    [X_trim, U_opt] = part_to_state_ctrl(gbest_pos, Va_ideal, psi_ideal);
end

% Graficas

function imprimir_resultados_pso(U_trim, J_trim, X_trim)
    Xd_trim = xdot(X_trim, U_trim);
    alpha_trim = atan2(X_trim(3), X_trim(1));
    Va_trim = norm(X_trim(1:3));
    gamma_trim = X_trim(8) - alpha_trim;
    sep = '=======================================================';

    fprintf('\n%s\n  OPTIMAL CONTROL SET\n%s\n', sep, sep);
    fprintf('  delta_A (aleron)    : %+12.6f deg  (%+14.6f rad)\n', rad2deg(U_trim(1)), U_trim(1));
    fprintf('  delta_E (elevador)  : %+12.6f deg  (%+14.6f rad)\n', rad2deg(U_trim(2)), U_trim(2));
    fprintf('  delta_R (timon)     : %+12.6f deg  (%+14.6f rad)\n', rad2deg(U_trim(3)), U_trim(3));
    fprintf('  throttle 1          : %29.8f\n', U_trim(4));
    fprintf('  throttle 2          : %29.8f\n', U_trim(5));
    fprintf('  Costo J final       : %29.6f\n', J_trim);

    fprintf('\n%s\n  TRIM STATE X_trim\n%s\n', sep, sep);
    labels = {'u', 'v', 'w', 'p', 'q', 'r', 'phi', 'theta', 'psi'};
    units = {'m/s', 'm/s', 'm/s', 'rad/s', 'rad/s', 'rad/s', 'rad', 'rad', 'rad'};
    for i = 1:9
        fprintf('  X[%d]  %-6s [%-5s] : %+18.8f\n', i-1, labels{i}, units{i}, X_trim(i));
    end

    fprintf('\n  -- Ángulos derivados --\n');
    fprintf('  alpha (ataque)      : %+12.6f deg  (%+14.6f rad)\n', rad2deg(alpha_trim), alpha_trim);
    fprintf('  gamma (senda)       : %+12.6f deg  (%+14.6f rad)\n', rad2deg(gamma_trim), gamma_trim);
    fprintf('  Va (aerodin.)       : %+12.6f m/s\n', Va_trim);

    fprintf('\n%s\n  XDOT evaluado en el trim\n%s\n', sep, sep);
    for i = 1:9
        fprintf('  Xd[%d] %-6s_dot      : %1.6e\n', i-1, labels{i}, Xd_trim(i));
    end
    fprintf('\n  ||Xdot|| = %1.6e  (ideal = 0)\n%s\n', norm(Xd_trim), sep);
end

function d = rad2deg(r), d = r * (180/pi); end
