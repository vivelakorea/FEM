% This file is implementation of Box2.5 in Belytschko. 

% 3 node 1D model

n_e = 2;
n_N = 3;

% mm, N 
rho = 0.0000027;
A0 = 100;
l(1,1) = 100; % l(element, time)
l(2,1) = 100;


time_step = 0.001;
% 1. Initial conditions and initialization: set v0, sig_e^0; n=1, t=time_step;
% modification: n starts at 1 
x(:,1) = [0;1;2];
u(:,1) = [0;0;0];
v(:,1) = [0;0;0]; % v(:,n) stands for v_(n+1/2)
a(:,1) = [0;0;0];

Le(:,:,1) = [1 0 0; 0 1 0];
Le(:,:,2) = [0 1 0; 0 0 1];

m1 = A0 * l(1,1) * rho / 6;
m2 = A0 * l(2,1) * rho / 6;
M = [2*m1 m1 0; m1 2*(m1+m2) m2; 0 m2 2*m2];
disp(inv(M))

n = 1;
t = time_step;

while (1)
    if n == 1000000
        break;
    end
    % 2. Get f(:,n)
    f(:, n) = [0;0;0];
    for i = 1:n_e % for every element
        % 1) Gather element nodal displacements u_e(n) and v_e(n+1/2)
        u_e(:,i) = Le(:,:,i) * u(:,n); % u_e([x1;x2], element)
        l_e = u_e(2,i) - u_e(1,i) + l(i,1);
        % 2) If n = 1, skip
        if n ~= 1
            % 3) Compute measure of deformation
            strain = (u_e(2,i) - u_e(1,i)) / l(i,1);
            fprintf("element: %d, strain: %f\n", i, strain)
            % 4) Compute stress by constitutive equation
            stress = swift(strain);
            fprintf("element: %d, stress: %f\n", i, stress)
        end
        % 5) Compute internal nodal forces by relevant equation
        A_e = A0 * (l(i,1) / l_e);
        if n == 1
            f_int_e = [0;0];
        else
            f_int_e = stress * A_e * [-1;1];
            % fprintf("element: %d, f_int_e: %f\n", i, f_int_e)
        end
            % 6) Compute external nodal forces on element and f_ext
        if i == 1
            f_ext_e = [-1;0];
        else
            f_ext_e = [0;1];
        end
        f_e = f_ext_e - f_int_e;
        % fprintf("element: %d, f_e: %f\n", i, f_e)
        % 7) Scatter element nodal forces to global matrices
        f(:,n) = f(:,n) + Le(:,:,i).' * f_e;
    end
    fprintf("time step: %d, u: [%f %f %f]'\n", n, u(1,n), u(1,n), u(2,n))
    fprintf("time step: %d, f: [%f %f %f]'\n", n, f(1,n), f(2,n), f(3,n))
    % 3. Compute accelartions a(:,n) = M\f(:,n)
    a(:,n) = M\f(:,n);
    fprintf("time step: %d, a: [%f %f %f]'\n", n, a(1,n), a(2,n), a(3,n))
    % 4. Update nodal velocities
    if n == 1
        v(:,1) = v(:,1) + 0.5 * time_step * a(:,1);
    else
        v(:,n) = v(:,n-1) + time_step * a(:,n);
    end
    % 5. Enforce essential boundary conditions
    v(1,n) = 0;
    fprintf("time step: %d, v: [%f %f %f]'\n", n, v(1,n), v(2,n), v(3,n))
    % 6. Update nodal displacements
    u(:,n+1) = u(:,n) + time_step * v(:,n);
    n = n + 1;
    t = t + time_step;

    xs = [u(1,n+1), u(2,n+1), u(3,n+1)];
    ys = [0,0,0];
    scatter(xs, ys)
    xlim([0 .1])
    pause(0.1)
end


function stress = swift(strain)
    stress = 2000 * (abs(strain + 1) ^ 0.2 - 1);
end