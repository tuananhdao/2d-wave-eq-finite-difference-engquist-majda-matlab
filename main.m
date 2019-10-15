clear
close all
clc

dt = 0.1;
dx = 2.5*dt;
dy = dx;
Nx = 51;
Ny = 121;

x = 0:dx:(Nx-1)*dx;y = 0:dy:(Ny-1)*dy;
[X,Y] = meshgrid(x,y);
r=10;r2=0.7;
w = exp(-(sqrt((X-x(end)).^2+(Y-y(end)).^2)-r).^2/r2)';
w_new = exp(-(sqrt((X-x(end)).^2+(Y-y(end)).^2)-r-dt).^2/r2)';

for n = 1 : 180
    w_old = w; w = w_new;
    
    pcolor(X,Y,w');title(n*dt);daspect([1 1 1]);colormap(hot);caxis([-.2,1]);colorbar;pause(0.01);
    
    % In the interior: w = 2w^n - w^(n-1) + (dt)^2(D2x w + D2y w)
    w_new(2:Nx-1,2:Ny-1) = 2*w(2:Nx-1,2:Ny-1)-w_old(2:Nx-1,2:Ny-1)...
                        +dt^2*( (w(3:Nx,2:Ny-1)-2*w(2:Nx-1,2:Ny-1)+w(1:Nx-2,2:Ny-1))/(dx^2)...
                                +(w(2:Nx-1,3:Ny)-2*w(2:Nx-1,2:Ny-1)+w(2:Nx-1,1:Ny-2))/(dy^2));

    w_new(end,:) = w_new(end-1,:); % right-hand wall
    w_new(2:Nx-1,end) = w_new(2:Nx-1,end-1); % upper wall
    w_new(2:Nx-1,1) = w_new(2:Nx-1,2); % lower wall
    
    % Perfectly reflecting Neumann
%     w_new(1,:) = w_new(2,:); % left-hand wall
%     w_new(1,:) = 0;
    % First approximation
%     w_new(1,:) = w(2,:)+(dx-dt)/(dx+dt)*w(1,:)+(-dx+dt)/(dx+dt)*w_new(2,:);
    % Second approximation
    w_new(1,:) = (2*dx*dt^2)/(dx+dt)*((1/2/dx/dt-1/2/dt^2-1/2/dy^2)*w_old(1,:)+(-1/2/dx/dt-1/2/dt^2)*w_old(2,:)...
        +(1/4/dy^2)*([w_old(1,2:end),w_old(1,end)]+[w_old(1,1),w_old(1,1:end-1)]...
        +[w_new(2,2:end),w_new(2,end)]+[w_new(2,1),w_new(2,1:end-1)])...
        +(1/dt^2)*w(1,:)+(1/dt^2)*w(2,:)...
        +(1/2/dx/dt-1/2/dt^2-1/2/dy^2)*w_new(2,:));
end