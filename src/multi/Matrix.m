function Matrix(m, L)
    syms a1 r11 r12 x1 da1 dr11 dr12 m1

    % Define position and velocity vectors
    R1 = [r11; r12];
    dR1 = [dr11; dr12];
    q1 = [R1; a1];
    dq1 = [dR1; da1];

    % Define kinetic energy
    T1 = (1/2) * m1 * (dR1.' * dR1);

    % Define rotation matrix A1
    A1 = [cos(a1), -sin(a1); sin(a1), cos(a1)];

    % Compute the derivative of the rotation matrix
    dA1 = diff(A1, a1);

    % Define position offset
    U_1 = [x1; 0];

    % Compute new position vector
    r1p = R1 + A1 * U_1;

    % Compute new velocity vector
    dr1p = dR1 + dA1 * U_1 * da1;

    % Define G1 matrix
    G1 = [eye(2), dA1 * U_1];

    % Display results
    disp('动能 T1:');
    disp(T1);

    disp('旋转矩阵 A1:');
    disp(A1);

    disp('旋转矩阵的导数 dA1:');
    disp(dA1);

    disp('新的位置矢量 r1p:');
    disp(r1p);

    disp('新的速度矢量 dr1p:');
    disp(dr1p);

    disp('G1 矩阵:');
    disp(G1);

    % Compute and simplify integrals
    M1_1 = simplify(m1 / L * int(G1.' * G1, x1, -0.5 * L, 0.5 * L));
    M1_2 = simplify(m1 / L * int(G1.' * G1, x1, -0.2 * L, 0.8 * L));

    % Display results
    disp('M1{1}:');
    disp(M1_1);

    disp('M1{2}:');
    disp(M1_2);
end
