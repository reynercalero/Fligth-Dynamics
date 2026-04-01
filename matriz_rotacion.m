function R = matriz_rotacion(phi, theta, psi)

    % Matriz de rotación roll (Eje X)
    Rx = [  1  ,    0    ,     0    ;
            0  , cos(phi), sin(phi);
            0  , -sin(phi), cos(phi)];

    % Matriz de rotación pitch (Eje Y)
    Ry = [cos(theta) ,  0 , -sin(theta) ;
               0     ,  1 ,     0      ;
          sin(theta),  0 , cos(theta)];

    % Matriz de rotación yaw (Eje Z)
    Rz = [cos(psi), sin(psi),  0   ;
          -sin(psi), cos(psi) ,  0   ;
             0    ,     0    ,  1  ];
    %matriz de rotacion en los 3 ejes
    R = Rz * Ry * Rx;
end
