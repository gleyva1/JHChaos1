function dot = threeBf(I,m1,m2,m3,G)
    % Authors: Gerard Leyva

    % Input and output.
    % Mass is entered through arguments m1,m2,m3. The G constant can be
    % set to 1 to make the computation unitless.
    %
    % The input I is expected to be a Matlab vector, with positions being
    % the first 6 elements (x1 position of mass 1, y1 position of mass 1,
    % x2 position of mass 2, etc.), followed by velocities for the
    % remaining six elements (vxi velocity of mass 1, vyi velocity of mass
    % 1, etc.).
    x1 = I(1);
    y1 = I(2);
    x2 = I(3);
    y2 = I(4);
    x3 = I(5);
    y3 = I(6);
    
    dot = zeros(1,12);

    % Physical parameters
    %G = 6.67384e-11;    % SI units
    %G = 0.498199;    % in m^3/(kg*day^2)
    
    %m2 = 5.97219e24;    % Mass of the Earth (in kg)
    %m3 = 7.34767309e22; % Mass of the Moon (in kg)
    %m1 = 1.9891e30;     % Mass of the Sun (in kg)

    % Separation distances rij
    r12 = sqrt((x1-x2)^2 + (y1-y2)^2);
    r23 = sqrt((x3-x2)^2 + (y3-y2)^2);
    r13 = sqrt((x1-x3)^2 + (y1-y3)^2);
    
    %{
    function Fij = gravF(mj,xi,xj,r)
        Fij = -G*mj*(xi-xj)/r^3;
    end
    %}

    %{
    % Dimensionalized form.
    function Fijk = grav(mj,mk,xi,xj,xk,rij,rik)
        Fijk = -G*mj*(xi-xj)/rij^3 - G*mk*(xi-xk)/rik^3;
    end
    %}

    % Once the dimensionless form is found, edit this part of the code.
    % The dimensionalized form is saved above.
    function Fijk = grav(mj,mk,xi,xj,xk,rij,rik)
        Fijk = -G*mj*(xi-xj)/rij^3 - G*mk*(xi-xk)/rik^3;
    end

    % ODEs. First positions, then velocities
    %{
    dot(1) = I(7);  % xdot1 = vx1
    dot(2) = I(8);  % ydot1 = vy1
    dot(3) = I(9);  % xdot2 = vx2
    dot(4) = I(10); % Ydot2 = vy2
    dot(5) = I(11); % xdot3 = vx3
    dot(6) = I(12); % ydot3 = vy3
    dot(7) = gravF(m2,x1,x2,r12) + gravF(m3,x1,x3,r13);   % vxdot1
    dot(8) = gravF(m2,y1,y2,r12) + gravF(m3,y1,y3,r13);   % vydot1
    dot(9) = gravF(m1,x2,x1,r12) + gravF(m3,x2,x3,r23);   % vxdot2
    dot(10) = gravF(m1,y2,y1,r12) + gravF(m3,y2,y3,r23);  % vydot2
    dot(11) = gravF(m1,x3,x1,r13) + gravF(m2,x3,x2,r23);  % vxdot3
    dot(12) = gravF(m1,y3,y1,r13) + gravF(m2,y3,y2,r23);  % vydot3
    %}

    dot(1) = I(7);  % xdot1 = vx1
    dot(2) = I(8);  % ydot1 = vy1
    dot(3) = I(9);  % xdot2 = vx2
    dot(4) = I(10); % Ydot2 = vy2
    dot(5) = I(11); % xdot3 = vx3
    dot(6) = I(12); % ydot3 = vy3
    dot(7) = grav(m2,m3,x1,x2,x3,r12,r13);   % vxdot1
    dot(8) = grav(m2,m3,y1,y2,y3,r12,r13);   % vydot1
    dot(9) = grav(m1,m3,x2,x1,x3,r12,r23);   % vxdot2
    dot(10) = grav(m1,m3,y2,y1,y3,r12,r23);  % vydot2
    dot(11) = grav(m1,m2,x3,x1,x2,r13,r23);  % vxdot3
    dot(12) = grav(m1,m2,y3,y1,y2,r13,r23);  % vydot3

end