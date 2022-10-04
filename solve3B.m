function S = solve3B(I,m1,m2,m3,G)
    % Authors: Gerard Leyva
    % To run this, put in the command window >>
    % I = [-0.5 0 0.5 0 -0.1 0.75 0 -0.3 0 0.3 0 -0.3];
    % solve3B(I,0.5,0.5,0.5,0.0075);
    t0 = 0;
    n = 2750000;
    hscale = 0.000002;

    % Solar system simulation. These scales measures time in days. The
    % other units are in usual SI. In threeBf.m, use lines 43-45 and
    % comment-out 50-56.
    % Type in console >>
    % SSday = [0 0 1.495979e11 0 1.495979e11 -3.84472e8 0 0 0 2.592e9 8.83008e7 2.592e9]
    % SEM = solve3B(SSday,1.9891e30,5.97219e24,7.34767309e22,0.498199);
    %G = 0.498199;    % in m^3/(kg*day^2);
    %m2 = 5.97219e24;    % Mass of the Earth (in kg)
    %m3 = 7.34767309e22; % Mass of the Moon (in kg)
    %m1 = 1.9891e30;     % Mass of the Sun (in kg)
    %n = 15000000;       % Uncomment this
    %hscale = 9^18;      % Uncomment this


    % Runge-Kutta function to solve the ODE
    %function R=rk4(t0,I,hscale,n,f)
    function R=rk4(t0,I,hscale,n)
        R=zeros(n+1,13);
        t=zeros(n+1,1);
        y=zeros(n+1,12);
        t(1)=t0;
        y(1,1)=I(1);      % x1
        y(1,2)=I(2);      % y1
        y(1,3)=I(3);      % x2
        y(1,4)=I(4);      % y2
        y(1,5)=I(5);      % x3
        y(1,6)=I(6);      % y3
        y(1,7)=I(7);      % vx1
        y(1,8)=I(8);      % vy2
        y(1,9)=I(9);      % vx2
        y(1,10)=I(10);    % vy2
        y(1,11)=I(11);    % vx3
        y(1,12)=I(12);    % vy3

        showpercentage = 1;
        for j=1:n
            r12 = sqrt((y(j,1)-y(j,3))^2 + (y(j,2)-y(j,4))^2);
            r23 = sqrt((y(j,5)-y(j,3))^2 + (y(j,6)-y(j,4))^2);
            r13 = sqrt((y(j,1)-y(j,5))^2 + (y(j,2)-y(j,6))^2);

            h = hscale/(r12^2 + r23^2 + r13^2);
            %h = hscale;

            %{
            %Original, 1D recipe
            t(j+1)=t(j)+h;
            k1=f(t(j),y(j));
            k2=f(t(j)+h/2,y(j)+h/2*k1);
            k3=f(t(j)+h/2,y(j)+h/2*k2);
            k4=f(t(j)+h,y(j)+h*k3);
            y(j+1)=y(j)+h/6*(k1+2*k2+2*k3+k4);
            %}

            % Modded for threeBf(I,m1,m2,m3,G)
            t(j+1)=t(j)+h;
            k1=threeBf(y(j,:),m1,m2,m3,G);
            k2=threeBf(y(j,:)+h/2*k1,m1,m2,m3,G);
            k3=threeBf(y(j,:)+h/2*k2,m1,m2,m3,G);
            k4=threeBf(y(j,:)+h*k3,m1,m2,m3,G);
            
            y(j+1,:)=y(j,:)+h/6*(k1+2*k2+2*k3+k4);

            if (100*j/n) > showpercentage
                disp([num2str(showpercentage),'% complete (t = ',num2str(t(j+1)),')'])
                showpercentage = showpercentage + 1;
            end

        end
        R(:,1)=t;
        R(:,2:13)=y;
    end

    S = rk4(t0,I,hscale,n);
    
    figure;close;hold on;
    plot(S(:,2),S(:,3),'r',S(:,4),S(:,5),'b',S(:,6),S(:,7),'k');
    %axis ([-3.5 3.5 -3.5 3.5]);
    scatter(S(1,2),S(1,3),'rd','filled');
    scatter(S(1,4),S(1,5),'bd','filled');
    scatter(S(1,6),S(1,7),'kd','filled');
    scatter(S(n+1,2),S(n+1,3),'r','filled');
    scatter(S(n+1,4),S(n+1,5),'b','filled');
    scatter(S(n+1,6),S(n+1,7),'k','filled');

    xlabel('$x(t)$','interpreter','latex');
    ylabel('$y(t)$','interpreter','latex');
    title({['Three Body System: $m_{1}$ in red, $m_{2}$ in blue, $m_{3}$ in black.'];['Starting positions: Diamonds | Final positions: Dots']},'interpreter','latex');

end