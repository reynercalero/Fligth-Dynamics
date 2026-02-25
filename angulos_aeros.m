function [alpha, beta, gamma] =angulos_aeros(Vb, theta)

    u = Vb(1); v = Vb(2); w = Vb(3);
    Vtotal = norm(Vb);

    % angulo de ataque (alpha)
    alpha =rad2deg(atan2(w, u));

    % angulo de sidelisp (beta)
    if Vtotal == 0
        beta = 0;
    else
        beta = rad2deg(asin(v / Vtotal));
    end

    %angulo de climb
    gamma = theta - alpha;
end
