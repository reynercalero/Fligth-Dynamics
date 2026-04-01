function visualizador(R, Vb,caso_vuelo)

    p_local = [1.5 ,   0  , 0;
                 0 ,  0.8 , 0;
                 0 , -0.8 , 0;
               -0.5,   0  , 0]';
    avion_rot = R * p_local;
    ejes_b = R * (eye(3) * 1); % Ejes del cuerpo rotados

    figure(1); clf;
    hold on; axis equal; grid off;

    view(0, 30);

    %ejes ned
    lim = 2.5;
    line([-lim, lim], [0, 0], [0, 0], 'Color', [0.6 0.6 0.6], 'LineStyle', '--'); text(lim,0,0,' N');
    line([0, 0], [-lim, lim], [0, 0], 'Color', [0.6 0.6 0.6], 'LineStyle', '--'); text(0,lim,0,' E');
    line([0, 0], [0, 0], [-lim, lim], 'Color', [0.6 0.6 0.6], 'LineStyle', '--'); text(0,0,lim,' D');

    %ejes body
    col_b = {'r', 'g', 'b'};
    for i = 1:3
        quiver3(0, 0, 0, ejes_b(1,i), ejes_b(2,i), ejes_b(3,i), ...
                'Color', col_b{i}, 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
    end

    %Dibujar el avion
    fill3(avion_rot(1,1:3), avion_rot(2,1:3), avion_rot(3,1:3), 'cyan', 'FaceAlpha', 0.3);
    plot3(avion_rot(1,[1,4]), avion_rot(2,[1,4]), avion_rot(3,[1,4]), 'k', 'LineWidth', 2);


    % vector de velocidad body
    Vb_rot = R * (Vb / (norm(Vb) + 1e-6) * 1.8);
    quiver3(0,0,0, Vb_rot(1), Vb_rot(2), Vb_rot(3), 'm', 'LineWidth', 1, 'MaxHeadSize', 0.2);
    text(Vb_rot(1), Vb_rot(2), Vb_rot(3), ' V_{body}', 'Color', 'm', 'FontWeight', 'bold');

    % ejes
    set(gca, 'zdir', 'reverse', 'ydir', 'reverse');
    title(['Gráfica de Transformación: Caso ', caso_vuelo]);
    xlabel('Norte'); ylabel('Este'); zlabel('Abajo');
    axis off;

end
