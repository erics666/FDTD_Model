% =========================================================================
% Zisheng Wang
% FDTD ABC condition 1D
% Final project of ECE 402
% =========================================================================

clc;clear;close all;

% =========================================================================
% CONSTANT
% =========================================================================
sigma_0 = 8.854187817e-12;
mu_0 = pi*4e-7;

% speed of light
c = 1/sqrt(sigma_0*mu_0);

% =========================================================================
% MAN SET PATAMETERS
% =========================================================================
% conductor of the material
sigma = 1 * sigma_0; % be air

% number of x
xNum = 400; 

% number of time step
tNum = 300;

% phase velocity
v_p = 1/sqrt(mu_0 * sigma);

% x step
delta_x = 5.4e-2; % example in textboot p.77

% time step
delta_t = 0.18e-9;

E_z = zeros(1, xNum); % plus 2 for ABC first order;
H_y = zeros(1, xNum);

E_z_record = zeros(tNum, xNum);
H_y_record = zeros(tNum, xNum);

% initialize with gaussian plus
E_in = exp( -((1:1:tNum)-10).^2/(10^2) );
plot((1:1:tNum)*delta_t, E_in);

E_mon = zeros(tNum,1);

h = figure(2);
filename = 'testAnimated.gif';

for tt = 1:1:tNum-1
    
    for xx = 1:1:xNum-1
        H_y(xx) = H_y(xx) + delta_t/(mu_0*delta_x) * ( E_z(xx+1) - E_z(xx) );
    end
    
    for xx = 1:1:xNum-1
        E_z(xx+1)   = E_z(xx+1) + delta_t/(sigma*delta_x) * ( H_y(xx+1) - H_y(xx) );
    end
        
    E_z(xNum/2) = E_in(tt);
    
    E_z_record(tt,:) = E_z;
    H_y_record(tt,:) = H_y;
    
    % Mur's Absorbing boundary condition
    if tt>1
        E_z(xNum)=dum2;%+((S-1)/(S+1))*(E_z(xNum-1)-E_z(xNum));
        E_z(2)=dum1;%+((S-1)/(S+1))*(E_z(2)-E_z(1));
        E_z(1)=E_z(2);
    end
    dum1=E_z(3);
    dum2=E_z(xNum-1);
    
%     figure(3);

    plot( (1:1:xNum)*delta_x, E_z);
    axis([0 xNum*delta_x -0.1 1])
    xlabel(['t=' num2str(tt*delta_t*10e6) '\mu s']);
    grid on;
    set(h,'Units','centimeters');
    set(h,'Position',[0 0 15 12*9/16]);
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if tt == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.03); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.03); 
    end 
    
    E_mon(tt) = abs( E_z(399)/abs(max(E_z)) ) * 100;
end

% h = figure(2);
% plot( (1:1:tNum)*delta_t*10e6, E_mon );
% xlabel('t [\mu s]'); ylabel('E_{mon} [%]')
% title('FDTD with absolute boundary conditions for 1D and TE problem');
% axis([ 0 tNum*delta_t*10e6 0 100 ]);
% grid on;
% set(h,'Units','centimeters');
% set(h,'Position',[0 0 15 12*9/16]);
