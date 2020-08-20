% =========================================================================
% Zisheng Wang
% FDTD ABC condition 2D
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
yNum = 400; 

% number of time step
tNum = 300;

xsource=yNum/2;
ysource=yNum/2;

% phase velocity
v_p = 1/sqrt(mu_0 * sigma);

S=1/(2^0.5);

% x step
delta_x = 5.4e-2; % example in textboot p.77

% time step
delta_t = S*delta_x/c;

E_z = zeros(yNum, yNum); % plus 2 for ABC first order;
H_x = zeros(yNum, yNum);
H_y = zeros(yNum, yNum);

E_z_record = zeros(tNum, yNum);
H_y_record = zeros(tNum, yNum);

% initialize with gaussian plus
E_in = exp( -((1:1:tNum)-10).^2/(10^2) );
plot((1:1:tNum)*delta_t, E_in);

E_mon = zeros(tNum,1);

p0=1;
p2=-0.5;
c0=(c/(2*S))*(1-(p0/S));
c1=-(c/(2*S))*(1+(p0/S));
c2=(c/(S^2))*(p0+(p2*S*S));
c3=-(p2*c)/2;
c0efffor=-(c0/c1);
c2efffor=-(c2/c1);
c3efffor=-(c3/c1);
c0=(c/(2*S))*(1+(p0/S));
c1=-(c/(2*S))*(1-(p0/S));
c2=-(c/(S^2))*(p0+(p2*S*S));
c3=(p2*c)/2;
c1effrev=-(c1/c0);
c2effrev=-(c2/c0);
c3effrev=-(c3/c0);

prev_xfor=zeros(1,yNum);
prev_x_minus_1for=zeros(1,yNum);
prev_yfor=zeros(yNum,1);
prev_y_minus_1for=zeros(yNum,1);
prev_xrev=zeros(1,yNum);
prev_x_minus_1rev=zeros(1,yNum);
prev_yrev=zeros(yNum,1);
prev_y_minus_1rev=zeros(yNum,1);

h = figure(2);
filename = 'testAnimated.gif';

for tt = 1:1:tNum-1
    
    % Setting time dependent boundaries to update only relevant parts of the 
    % vector where the wave has reached to avoid unnecessary updates.
    
    n1 = 2; n2 = yNum-2; n11 = 2; n22 = yNum-2;

    %Vector update instead of for-loop for Hy and Hx fields
    H_x(n1:n2-1,n11:n22-1)=H_x(n1:n2-1,n11:n22-1)-delta_t/(mu_0*delta_x)*(E_z(n1:n2-1,n11+1:n22)-E_z(n11:n2-1,n11:n22-1));
    H_y(n1:n2-1,n11:n22-1)=H_y(n1:n2-1,n11:n22-1)+delta_t/(mu_0*delta_x)*(E_z(n1+1:n2,n11:n22-1)-E_z(n11:n2-1,n11:n22-1));
    
    %Vector update instead of for-loop for Ez field
    E_z(n1+1:n2-1,n11+1:n22-1)=E_z(n1+1:n2-1,n11+1:n22-1)+(H_y(n1+1:n2-1,n11+1:n22-1)-H_y(n1:n2-2,n11+1:n22-1)-H_x(n1+1:n2-1,n11+1:n22-1)+H_x(n1+1:n2-1,n11:n22-2))*delta_t/(sigma*delta_x);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %forward boundary
    if tt>=yNum-2-xsource
        E_z(yNum-2,3:1:yNum-3)=c0efffor*(E_z(yNum-3,3:1:yNum-3)+prev_prev_xfor(1,3:1:yNum-3))-prev_prev_x_minus_1for(1,3:1:yNum-3)+c2efffor*(prev_xfor(1,3:1:yNum-3)+prev_x_minus_1for(1,3:1:yNum-3))+c3efffor*(prev_x_minus_1for(1,2:1:yNum-4)+prev_x_minus_1for(1,4:1:yNum-2)+prev_xfor(1,2:1:yNum-4)+prev_xfor(1,4:1:yNum-2));
    end
    
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at forward boundary
    prev_prev_xfor=prev_xfor;
    prev_prev_x_minus_1for=prev_x_minus_1for;
    prev_xfor(1,1:1:yNum)=E_z(yNum-2,1:1:yNum);
    prev_x_minus_1for(1,1:1:yNum)=E_z(yNum-3,1:1:yNum);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %backward boundary
    if tt>=xsource-3
        E_z(2,3:1:yNum-3)=-prev_prev_xrev(1,3:1:yNum-3)+c1effrev*(E_z(3,3:1:yNum-3)+prev_prev_x_minus_1rev(1,3:1:yNum-3))+c2effrev*(prev_xrev(1,3:1:yNum-3)+prev_x_minus_1rev(1,3:1:yNum-3))+c3effrev*(prev_x_minus_1rev(1,2:1:yNum-4)+prev_x_minus_1rev(1,4:1:yNum-2)+prev_xrev(1,2:1:yNum-4)+prev_xrev(1,4:1:yNum-2));
    end
    
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at backward boundary
    prev_prev_xrev=prev_xrev;
    prev_prev_x_minus_1rev=prev_x_minus_1rev;
    prev_xrev(1,1:1:yNum)=E_z(3,1:1:yNum);
    prev_x_minus_1rev(1,1:1:yNum)=E_z(2,1:1:yNum);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %upward boundary
    if tt>=yNum-2-ysource
        E_z(3:1:yNum-3,yNum-2)=c0efffor*(E_z(3:1:yNum-3,yNum-3)+prev_prev_yfor(3:1:yNum-3,1))-prev_prev_y_minus_1for(3:1:yNum-3,1)+c2efffor*(prev_yfor(3:1:yNum-3,1)+prev_y_minus_1for(3:1:yNum-3,1))+c3efffor*(prev_y_minus_1for(2:1:yNum-4,1)+prev_y_minus_1for(4:1:yNum-2,1)+prev_yfor(2:1:yNum-4,1)+prev_yfor(4:1:yNum-2,1));
    end
    
    yNum = 1;
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at upward boundary
    prev_prev_yfor=prev_yfor;
    prev_prev_y_minus_1for=prev_y_minus_1for;
    prev_yfor(1:1:yNum,1)=E_z(1:1:yNum,yNum-2);
    prev_y_minus_1for(1:1:yNum,1)=E_z(1:1:yNum,yNum-3);
    
    %Mur's abc conditions obtained from Mur's difference equation for
    %downward boundary
    if tt>=ysource-3
        E_z(3:1:yNum-3,2)=-prev_prev_yrev(3:1:yNum-3,1)+c1effrev*(E_z(3:1:yNum-3,3)+prev_prev_y_minus_1rev(3:1:yNum-3,1))+c2effrev*(prev_yrev(3:1:yNum-3,1)+prev_y_minus_1rev(3:1:yNum-3,1))+c3effrev*(prev_y_minus_1rev(2:1:yNum-4,1)+prev_y_minus_1rev(4:1:yNum-2,1)+prev_yrev(2:1:yNum-4,1)+prev_yrev(4:1:yNum-2,1));
    end
       
    %Storage vectors for boundary and boundary-1 values of previous and its
    %previous time steps updated at downward boundary
    prev_prev_yrev=prev_yrev;
    prev_prev_y_minus_1rev=prev_y_minus_1rev;
    prev_yrev(1:1:yNum,1)=E_z(1:1:yNum,3);
    prev_y_minus_1rev(1:1:yNum,1)=E_z(1:1:yNum,2);
        
    %Mirroring of corner values taking the fact that corners are reached by the fields from the previous corners
    %in two time steps as S=1/sqrt(2) viz. sqrt(2)*delta(distance between two corners) is reached in 2 time steps
    E_z(2,2)=prev_prev_xrev(3);
    E_z(2,yNum-2)=prev_prev_xrev(yNum-3);
    E_z(yNum-2,2)=prev_prev_x_minus_1for(3);
    E_z(yNum-2,yNum-2)=prev_prev_x_minus_1for(yNum-3);
    
    % Source conditions
    E_z(xsource,ysource)=E_in(tt);

    
    %Movie type colour scaled image plot of Ez
    imagesc(delta*(1:1:yNum)*1e+6,(1e+6*delta*(1:1:yNum))',E_z',[-1,1]);colorbar;
    title(['\fontsize{20}Colour-scaled image plot of Ez in a spatial domain with Mur absorbing boundary and at time = ',num2str(round(tt*deltat*1e+15)),' fs']); 
    xlabel('x (in um)','FontSize',20);
    ylabel('y (in um)','FontSize',20);
    set(gca,'FontSize',20);
    getframe;
   
%     set(h,'Position',[0 0 15 12*9/16]);
    % Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if tt == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.03); 
%     else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.03); 
%     end 
%     
%     E_mon(tt) = abs( E_z(399)/abs(max(E_z)) ) * 100;
end

% h = figure(2);
% plot( (1:1:tNum)*delta_t*10e6, E_mon );
% xlabel('t [\mu s]'); ylabel('E_{mon} [%]')
% title('FDTD with absolute boundary conditions for 1D and TE problem');
% axis([ 0 tNum*delta_t*10e6 0 100 ]);
% grid on;
% set(h,'Units','centimeters');
% set(h,'Position',[0 0 15 12*9/16]);
