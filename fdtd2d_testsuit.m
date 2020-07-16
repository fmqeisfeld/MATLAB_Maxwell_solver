function main
clear all;
clc;

global xdim ydim
global material 
global Ezx Ezy Hy Hx Jz Pinst 
global Ex Ey Hzx Hzy Jx Jy
global Jzx Jzy Jlzx Jlzy Jlx Jly
global Pzx Pzy Px Py
global time 

%spatially interpolated 'old' current density fields
%needed for the calculation of instantaneously absorbed
%power density
global Jx_prime Jy_prime Jz_prime Jlz_prime Jlx_prime Jly_prime 

global ezx1 ezy1 ezx2 ezy2 ez3  % coefficient matrix for E-field update
global hx1 hy1 hx2 hy2

global epsilon mu eps0 mu0 c sigma_stary sigma_starx sigmax sigmay
global bound_width
global x1 x2 y1 y2
global dt delta
global Imp0
global epsinf omegap gammap omegapL omega0L gammaL

global frequency N_lambda2d
global E0 %peak e-field strenght (from peak intensity)
global T_e T_i 
global nevery
global t0 sigma_t srcx srcy srcw src_spatial %laser params   

global n1 n2 n11 n21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Drude-Lorentz parameters for 'cold' Aluminum
epsinf  = 2.73;
omegap  = 2.2955e16;
gammap  = 1.1174e15;
omegapL = 7.6595e15;
omega0L = 2.4024e15;
gammaL  = 4.5199e14;
%%%%%%%%%%%%%%%%%%%%% NATURAL CONSTANTS  %%%%%%%%%%%%%%%%%%%%%%%%%
eps0=(1/(36*pi))*1e-9;  %Farad/m
mu0=4*pi*1e-7;  % Henrys per meter
c=2.998e8;        % m/s
Imp0 = sqrt(mu0 / eps0); % Vacuum Impedance

%%%%%%%%%%%%%%%% SIMLATION AND SAMPLE DIMENSIONS  %%%%%%%%%%%%%%%%%%%%%%%
xdim=600;   %Simulation box dimensions
ydim=600;

time=0;      %starting simulation time in seconds
tmax=4e-13;  %total simulation time in s
nevery=10;   %how often to plot? (every nevery integration steps)

bound_width=15; % perfectly matched layers thickness in units of FD-cells
delta=10e-9     % spatial discretization of FD-grid (x and y direction) in meters

%%%%%%%%%%%%%%%%%%%%%% SAMPLE GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%
n1=1;       %aux. varis for vectorized update eqns.
n2=xdim-1;
n11=1;    
n21=ydim-1;

x1=350;     %start x-coordinate of rectangular sample
x2=x1+150;  %end pos ...

y1=bound_width+20;      %...
y2=ydim-bound_width-20; %:..

material=zeros(xdim,ydim); 
T_e=zeros(ydim,xdim);
T_i=zeros(ydim,xdim);

for i=1:xdim
    for j=1:ydim
        %rectangular bulk
        if i>=x1 && i<=x2 && j>=y1 && j<=y2
            material(i,j)=1;            
        end                   
        %cavity
        if sqrt((i-x1).^2 + (j-y2/2).^2) < 100
            material(i,j)=0;
        end
        %droplets
        if sqrt((i-200).^2 + (j-y2/2).^2) < 30
            material(i,j)=1;
        end                
        if sqrt((i-300).^2 + (j-y2/2).^2) < 30
            material(i,j)=1;
        end
        if sqrt((i-250).^2 + (j-y2/2+80).^2) < 20
            material(i,j)=1;
        end        
        if sqrt((i-250).^2 + (j-y2/2-80).^2) < 20
            material(i,j)=1;
        end                
        
    end
end 

T_e=material.*300;
T_i=material.*300;

%check geometry
%imagesc((1:1:xdim),((1:1:ydim))',T_e');colorbar;
%return

%%%%%%%%%%%%%%%%%%%  LASER PARMETERS  %%%%%%%%%%%%%%%%%
I0=5e15;             % Incident laser intensity in W/m^2
lambda=800e-9;      %central wavelength im meters

srcx=50;       %x-coordinate of laser source in units of FD-cells
srcw=80;       %beam radius in units of FD-cells

srcy=round(ydim/2);% y-coordinate of laser source in units of FD-cells
%Spatial shape of laser beam (gaussian here)
src_spatial=exp(-0.5.*([bound_width:ydim-bound_width]-srcy).^2./srcw.^2);  
t0=50e-15;     %time of peak intensity of gaussian pulse in seconds
t_FWHM=50e-15; %FWHM of pulse duration in s

%----------------------------------------------------
E0=sqrt(2.0*I0*Imp0) %peak efield
E0=E0*2.0;           %because wave splits into two waves
E0=E0/sqrt(6);       %modulus of electric field strength is evenly 
                     %distributed on TMZ and TEZ modes
frequency = c/lambda;
N_lambda2d=c/(frequency*delta) %points per lambda. Should be > 20 at least
sigma_t=t_FWHM/2.355;  %sigma of pulse duration for gaussian pulse

%%%%%%%%%%%%%%%%%%% SOLVER PARAMS %%%%%%%%%%%%%%%%%%%
S=1/(2^0.5);        %Courant stability factor
S=0.95*S;           %
dt=S*delta/c        %integration time-step
nmax=round(tmax/dt) %max nr of steps  

%%%%%%%%%%%%%%%%%%  PREPARE FIELD VARIABLES %%%%%%%%%%%%%%%%%%%%
Ezx=zeros(xdim,ydim);
Ezy=zeros(xdim,ydim);
Hy=zeros(xdim,ydim);
Hx=zeros(xdim,ydim);

Ex=zeros(xdim,ydim);
Ey=zeros(xdim,ydim);
Hzx=zeros(xdim,ydim);
Hzy=zeros(xdim,ydim);
Jx=zeros(xdim,ydim);
Jy=zeros(xdim,ydim);

Jz=zeros(xdim,ydim); 
Jzx=zeros(xdim,ydim);
Jzy=zeros(xdim,ydim);

Jlzx=zeros(xdim,ydim);
Jlzy=zeros(xdim,ydim);
Jlx=zeros(xdim,ydim);
Jly=zeros(xdim,ydim);

Pzx=zeros(xdim,ydim);
Pzy=zeros(xdim,ydim);
Px=zeros(xdim,ydim);
Py=zeros(xdim,ydim);

% 'old' current densities (interpolated) 
%  needed to compute power density
Jx_prime  =zeros(xdim,ydim);
Jy_prime  =zeros(xdim,ydim);
Jz_prime  =zeros(xdim,ydim);
Jlz_prime =zeros(xdim,ydim);
Jlx_prime =zeros(xdim,ydim);
Jly_prime =zeros(xdim,ydim);

Pinst=zeros(xdim,ydim); %instantenously absorbed power density

%prepare coeff matrices for update equations
ezx1=zeros(xdim,ydim);
ezx2=zeros(xdim,ydim);
ezy1=zeros(xdim,ydim);
ezy2=zeros(xdim,ydim);
ez3=zeros(xdim,ydim);
hx1=zeros(xdim,ydim);
hx2=zeros(xdim,ydim);
hy1=zeros(xdim,ydim);
hy2=zeros(xdim,ydim);

% Initializing electric conductivity matrices in x and y directions
sigmax=zeros(xdim,ydim);      % electrical conductivity in Siemens/meters
sigmay=zeros(xdim,ydim);

sigma_starx=zeros(xdim,ydim); % 'magnetic' conductivity 
sigma_stary=zeros(xdim,ydim);

%%%%%%%%%%%%%%%%%%%%%%%%% COEFFS FOR PML  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradingorder=6;     %Order of polynomial on which sigma is modeled
refl_coeff=1e-6;    %Required reflection co-efficient

sigmamax=(-log10(refl_coeff)*(gradingorder+1)*eps0*c)/(2*bound_width*delta);
bf = sigmamax / (bound_width^gradingorder * (gradingorder+1));

x=0;
for j=bound_width+1:-1:1                
    for i=1:xdim
        sigmax(i,j)=         bf*((x+0.5)^(gradingorder+1)-(x-0.5*(j<bound_width+1))^(gradingorder+1));   
        sigmax(i,ydim-(j-1))=bf*((x+0.5)^(gradingorder+1)-(x-0.5*(j<bound_width+1))^(gradingorder+1));
    end    
    x=x+1;
end

x=0;
for i=bound_width+1:-1:1
    for j=1:ydim
        sigmay(i,j)=         bf*((x+0.5)^(gradingorder+1)-(x-0.5*(x>0))^(gradingorder+1));   
        sigmay(xdim-(i-1),j)=bf*((x+0.5)^(gradingorder+1)-(x-0.5*(x>0))^(gradingorder+1));
    end
    x=x+1;    
end

%Magnetic conductivity matrix obtained by Perfectly Matched Layer condition
%This is also split into x and y directions in Berenger's model
for i=1:xdim
    for j=1:ydim
        sigma_starx(i,j)=sigmax(i,j)*mu0/eps0;
        sigma_stary(i,j)=sigmay(i,j)*mu0/eps0;
    end
end

for i=1:xdim
    for j=1:ydim
        if material(i,j) > 0
            ezx1(i,j)=1.0;
            ezy1(i,j)=1.0;
            
            ezx2(i,j)=dt/eps0/epsinf/delta;
            ezy2(i,j)=dt/eps0/epsinf/delta;
            
            ez3(i,j)=dt/eps0/epsinf;
            
        else
            ezx1(i,j)=(eps0-0.5*dt*sigmax(i,j))./(eps0+0.5*dt*sigmax(i,j));
            ezy1(i,j)=(eps0-0.5*dt*sigmay(i,j))./(eps0+0.5*dt*sigmay(i,j));
            
            ezx2(i,j)=dt./(eps0+0.5*dt*sigmax(i,j))./delta;
            ezy2(i,j)=dt./(eps0+0.5*dt*sigmay(i,j))./delta;
                                         
        end
        
        hx1(i,j)=(mu0-0.5*dt*sigma_starx(i,j))/(mu0+0.5*dt*sigma_starx(i,j));
        hy1(i,j)=(mu0-0.5*dt*sigma_stary(i,j))/(mu0+0.5*dt*sigma_stary(i,j));
        
        hx2(i,j)=dt/(mu0+0.5*dt*sigma_starx(i,j))/delta;
        hy2(i,j)=dt/(mu0+0.5*dt*sigma_stary(i,j))/delta;
        
    end
end

hfig=figure('Name','FDTD-2D demo','units','normalized','outerposition',[0 0 1 1]);
axis tight manual
imgcnt=0; %image counter for filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:1:nmax
    do_Maxwell(n);            
    do_Diffusion();
    time=time+dt;
    
    if mod(n,nevery)==0
        
        % PLOTS
        subplot(2,2,1);        
        %E-field magnitude
        imagesc(delta*(1:1:xdim),(delta*(1:1:ydim))',(sqrt(Ezx.^2+Ezy.^2+Ex.^2+Ey.^2))');
        c1=colorbar;
        
        caxis([-1 1].*E0);
        xlabel('x in m');
        ylabel('y in m');
        c1.Label.String='|E| in V/m';
        
        %txt=strcat('t= ', num2str(time));
        txt={strcat('t:', num2str(time)), strcat('n:', num2str(n))};
        text(delta*xdim*0.75, delta*ydim*0.8,txt);
        
        %power density
        subplot(2,2,2)
        imagesc(delta*(1:1:xdim),(delta*(1:1:ydim))',(Pinst)');
        xlabel('x in m');
        ylabel('y in m');
        c4=colorbar;
        c4.Label.String='Power density in W/m^3';        
        
        %Electron temperature distribution
        subplot(2,2,3)
        imagesc(delta*(1:1:xdim),(delta*(1:1:ydim))',(T_e)');
        xlabel('x in m');
        ylabel('y in m');
        c2=colorbar;
        c2.Label.String='T_e in K';
        
        %Lattice temperature
        subplot(2,2,4)
        imagesc(delta*(1:1:xdim),(delta*(1:1:ydim))',(T_i)');
        xlabel('x in m');
        ylabel('y in m');
        c3=colorbar;
        c3.Label.String='T_i in K';
                
        frame=getframe(hfig);                        
        
%         %write frames to files ?
%         im = frame2im(frame);
%         [imind,cm] = rgb2ind(im,256);
%         
%         
%         if n==1 || mod(n,250)==0
%             imgcnt = imgcnt+1;
%             fname=strcat('my',num2str(imgcnt),'.jpg');
%             imwrite(imind,cm,fname,'jpg');
%         end
        
      
    end

end

end


                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO MAXWELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_Maxwell(n)
    global eps0 mu0 
    global Hy Hx Ezx Ezy Ez  
    global Ex Ey Jx Jy Hzx Hzy
    
    global xdim ydim
    global dt delta                
    global Jzx Jzy 
    global Jlzx Jlzy Jlx Jly
    global Pzx Pzy Px Py
   
    global n1 n2 n11 n21 
    global ezx1 ezy1 ezx2 ezy2 ez3
    global jz1 jz2 
    global hx1 hx2 hy1 hy2

    global gammaL omegapL omega0L gammap omegap
    global material
    
    global Jx_prime Jy_prime Jz_prime Jlz_prime Jlx_prime Jly_prime 
    global Pinst
    
    do_softsource2d(n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       %Ezx & Ezy update    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ezx(n1+1:n2+1,n11+1:n21+1)=ezx1(n1+1:n2+1,n11+1:n21+1).*Ezx(n1+1:n2+1,n11+1:n21+1) + ...
                               ezx2(n1+1:n2+1,n11+1:n21+1).*(Hy(n1+1:n2+1,n11+1:n21+1)-Hy(n1:n2,n11+1:n21+1)) - ...
                                ez3(n1+1:n2+1,n11+1:n21+1).*(Jzx(n1+1:n2+1,n11+1:n21+1)+Jlzx(n1+1:n2+1,n11+1:n21+1));

    Ezy(n1+1:n2+1,n11+1:n21+1)=ezy1(n1+1:n2+1,n11+1:n21+1).*Ezy(n1+1:n2+1,n11+1:n21+1) - ...
                               ezy2(n1+1:n2+1,n11+1:n21+1).*(Hx(n1+1:n2+1,n11+1:n21+1)-Hx(n1+1:n2+1,n11:n21)) - ...
                                ez3(n1+1:n2+1,n11+1:n21+1).*(Jzy(n1+1:n2+1,n11+1:n21+1)+Jlzy(n1+1:n2+1,n11+1:n21+1));                               

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       %Ex & Ey update      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hz=Hzx+Hzy;
    %ex1 = ezy1, ex2= ezy2
    %ey1 = ezx1, ey2= etx2
    Ex(n1+1:n2+1,n11+1:n21+1)=ezy1(n1+1:n2+1,n11+1:n21+1).*Ex(n1+1:n2+1,n11+1:n21+1) + ...
                              ezy2(n1+1:n2+1,n11+1:n21+1).*(Hz(n1+1:n2+1,n11+1:n21+1)-Hz(n1+1:n2+1,n11:n21)) - ...
                               ez3(n1+1:n2+1,n11+1:n21+1).*(Jx(n1+1:n2+1,n11+1:n21+1)+Jlx(n1+1:n2+1,n11+1:n21+1));
                           
    Ey(n1+1:n2+1,n11+1:n21+1)=ezx1(n1+1:n2+1,n11+1:n21+1).*Ey(n1+1:n2+1,n11+1:n21+1) - ...
                              ezx2(n1+1:n2+1,n11+1:n21+1).*(Hz(n1+1:n2+1,n11+1:n21+1)-Hz(n1:n2,n11+1:n21+1)) - ...
                               ez3(n1+1:n2+1,n11+1:n21+1).*(Jy(n1+1:n2+1,n11+1:n21+1)+Jly(n1+1:n2+1,n11+1:n21+1));    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Pzx,Pzy,Px,Py update    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pzx(n1+1:n2+1,n11+1:n21+1)=Pzx(n1+1:n2+1,n11+1:n21+1) + dt.* Jlzx(n1+1:n2+1,n11+1:n21+1);
    Pzy(n1+1:n2+1,n11+1:n21+1)=Pzy(n1+1:n2+1,n11+1:n21+1) + dt.* Jlzy(n1+1:n2+1,n11+1:n21+1);
    Px(n1+1:n2+1,n11+1:n21+1)=Px(n1+1:n2+1,n11+1:n21+1) + dt.* Jlx(n1+1:n2+1,n11+1:n21+1);
    Py(n1+1:n2+1,n11+1:n21+1)=Py(n1+1:n2+1,n11+1:n21+1) + dt.* Jly(n1+1:n2+1,n11+1:n21+1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Update Jz & Jlz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    j1=(2-dt*gammap)/(2+dt*gammap);
    j2=2*dt/(2+dt*gammap)*eps0*omegap^2;   
            
    Jzx(n1+1:n2+1,n11+1:n21+1)=j1.*Jzx(n1+1:n2+1,n11+1:n21+1) + ...
                              material(n1+1:n2+1,n11+1:n21+1).* ...
                              j2.*Ezx(n1+1:n2+1,n11+1:n21+1);
                          
    Jzy(n1+1:n2+1,n11+1:n21+1)=j1.*Jzy(n1+1:n2+1,n11+1:n21+1) + ...
                               material(n1+1:n2+1,n11+1:n21+1).* ...
                               j2.*Ezy(n1+1:n2+1,n11+1:n21+1);
                          
    %Update coeffs
    jl1=(2-dt*gammaL)/(2+dt*gammaL);
    jl2=2*dt/(2+dt*gammaL).*(omegapL^2*eps0);
    jl3=2*dt/(2+dt*gammaL).*(omega0L^2);
    Jlzx(n1+1:n2+1,n11+1:n21+1)=jl1.*Jlzx(n1+1:n2+1,n11+1:n21+1) + ...
                                material(n1+1:n2+1,n11+1:n21+1).* ...
                                jl2.*Ezx(n1+1:n2+1,n11+1:n21+1) - ...
                                jl3.*Pzx(n1+1:n2+1,n11+1:n21+1);

    Jlzy(n1+1:n2+1,n11+1:n21+1)=jl1.*Jlzy(n1+1:n2+1,n11+1:n21+1) + ...
                                material(n1+1:n2+1,n11+1:n21+1).* ...        
                                jl2.*Ezy(n1+1:n2+1,n11+1:n21+1) - ...
                                jl3.*Pzy(n1+1:n2+1,n11+1:n21+1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Update Jx, Jy & Jlx, Jly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Jx(n1+1:n2+1,n11+1:n21+1)=j1.*Jx(n1+1:n2+1,n11+1:n21+1) + ...
                              material(n1+1:n2+1,n11+1:n21+1).* ...        
                              j2.*Ex(n1+1:n2+1,n11+1:n21+1);

    Jy(n1+1:n2+1,n11+1:n21+1)=j1.*Jy(n1+1:n2+1,n11+1:n21+1) + ...
                              material(n1+1:n2+1,n11+1:n21+1)>0.* ...                
                              j2.*Ey(n1+1:n2+1,n11+1:n21+1);

    Jlx(n1+1:n2+1,n11+1:n21+1)=jl1.*Jlx(n1+1:n2+1,n11+1:n21+1) + ... 
                               material(n1+1:n2+1,n11+1:n21+1).* ...                
                               jl2.*Ex(n1+1:n2+1,n11+1:n21+1) - ...
                               jl3.*Px(n1+1:n2+1,n11+1:n21+1);   
                           
    Jly(n1+1:n2+1,n11+1:n21+1)=jl1.*Jly(n1+1:n2+1,n11+1:n21+1) + ...
                               material(n1+1:n2+1,n11+1:n21+1).* ...                
                               jl2.*Ey(n1+1:n2+1,n11+1:n21+1) - ...
                               jl3.*Py(n1+1:n2+1,n11+1:n21+1);                              
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Hx & Hy update    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    Ez=Ezx+Ezy;
    Hx(n1:n2,n11:n21)=hx1(n1:n2,n11:n21).*Hx(n1:n2,n11:n21)- ...
                      hx2(n1:n2,n11:n21).*(Ez(n1:n2,n11+1:n21+1)-Ez(n1:n2,n11:n21));
        
    Hy(n1:n2,n11:n21)=hy1(n1:n2,n11:n21).*Hy(n1:n2,n11:n21)+ ...
                      hy2(n1:n2,n11:n21).*(Ez(n1+1:n2+1,n11:n21)-Ez(n1:n2,n11:n21));              

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Hzx & Hzy update    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
    %hzx1 == hy1, hzx2=hy2, hzy1=hx1, hzy2=hx2;    
    Hzx(n1:n2,n11:n21)=hy1(n1:n2,n11:n21).*Hzx(n1:n2,n11:n21)- ...
                       hy2(n1:n2,n11:n21).*(Ey(n1+1:n2+1,n11:n21)-Ey(n1:n2,n11:n21));
                   
    Hzy(n1:n2,n11:n21)=hx1(n1:n2,n11:n21).*Hzy(n1:n2,n11:n21)+ ...
                       hx2(n1:n2,n11:n21).*(Ex(n1:n2,n11+1:n21+1)-Ex(n1:n2,n11:n21));
                   
    %Eabs=sqrt(Ezx.^2+Ezy.^2+Ex.^2+Ey.^2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute power density   
    % Current densities need to be interpolated (spatially and temporal)
    % to be consistent with E-field 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Jx_intp = 0.5.*((Jx(n1+1:n2+1, n11+1:n21+1)  + Jx(n1+1:n2+1,n11:n21)).*0.5 + Jx_prime(n1+1:n2+1,n11+1:n21+1));
    Jy_intp = 0.5.*((Jy(n1+1:n2+1, n11+1:n21+1)  + Jy(n1:n2,n11+1:n21+1)).*0.5 + Jy_prime(n1+1:n2+1,n11+1:n21+1));
    Jz_intp = 0.5.*((Jzx(n1+1:n2+1, n11+1:n21+1)  + Jzy(n1+1:n2+1,n11+1:n21+1)).*0.5 + Jz_prime(n1+1:n2+1,n11+1:n21+1));
    
    Jlx_intp = 0.5.*((Jlx(n1+1:n2+1, n11+1:n21+1)  + Jlx(n1+1:n2+1,n11:n21)).*0.5 + Jlx_prime(n1+1:n2+1,n11+1:n21+1));
    Jly_intp = 0.5.*((Jly(n1+1:n2+1, n11+1:n21+1)  + Jly(n1:n2,n11+1:n21+1)).*0.5 + Jly_prime(n1+1:n2+1,n11+1:n21+1));
    Jlz_intp = 0.5.*((Jlzx(n1+1:n2+1, n11+1:n21+1)  + Jlzy(n1+1:n2+1,n11+1:n21+1)).*0.5 + Jlz_prime(n1+1:n2+1,n11+1:n21+1));
    
    %prefactors
    foo = gammap./eps0./omegap.^2;
    bar = gammaL./eps0./omegapL.^2;
    %TMZ
    q_Drude   = foo .* Jz_intp.*Jz_intp;
    q_Lorentz = bar .* Jlz_intp.*Jlz_intp;
    %TEZ
    q_Drude   = q_Drude   + foo .* (Jx_intp.*Jx_intp + Jy_intp.*Jy_intp);
    q_Lorentz = q_Lorentz + bar .* (Jlx_intp.*Jlx_intp + Jly_intp.*Jly_intp);
    
    Pinst(n1:n2,n11:n21) = q_Drude + q_Lorentz;
    
    %save 'old' current densities for next iteration step
    %Drude
    Jx_prime(n1+1:n2+1, n11+1:n21+1) = (Jx(n1+1:n2+1, n11+1:n21+1)  + Jx(n1+1:n2+1,n11:n21)).*0.5;
    Jy_prime(n1+1:n2+1, n11+1:n21+1) = (Jy(n1+1:n2+1, n11+1:n21+1)  + Jy(n1:n2,n11+1:n21+1)).*0.5;
    Jz_prime(n1+1:n2+1, n11+1:n21+1) = (Jzx(n1+1:n2+1, n11+1:n21+1)  + Jzy(n1+1:n2+1,n11+1:n21+1)).*0.5;
    %Lorentz
    Jlx_prime(n1+1:n2+1, n11+1:n21+1) = (Jlx(n1+1:n2+1, n11+1:n21+1)  + Jlx(n1+1:n2+1,n11:n21)).*0.5;
    Jly_prime(n1+1:n2+1, n11+1:n21+1) = (Jly(n1+1:n2+1, n11+1:n21+1)  + Jly(n1:n2,n11+1:n21+1)).*0.5;
    Jlz_prime(n1+1:n2+1, n11+1:n21+1) = (Jlzx(n1+1:n2+1, n11+1:n21+1)  + Jlzy(n1+1:n2+1,n11+1:n21+1)).*0.5;
    
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% softSource 2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_softsource2d(n)
    global Ezx Ezy Hzx Hzy 
    global delta
    global dt
    global c eps0 mu0 Imp0
    global E0 %peak e-field strength    
    global N_lambda2d
    global ydim     
    global bound_width
    global t0 sigma_t srcx srcy src_spatial
            
    timefun=2*E0*exp(-0.5*(n*dt-t0)^2/sigma_t^2);        % Gauss
    sine=sin(((2*pi*(c/(delta*N_lambda2d))*n*dt))); 
                           
    Einc=timefun.*sine.*src_spatial;    %incident E-field amp
    Hinc=Einc./Imp0;                    %incident H-field amp
    %TMZ
    Ezx(srcx,bound_width:ydim-bound_width)=Ezx(srcx,bound_width:ydim-bound_width)+dt./eps0./delta.*Hinc;
    Ezy(srcx,bound_width:ydim-bound_width)=Ezy(srcx,bound_width:ydim-bound_width)+dt./eps0./delta.*Hinc;
    %TEZ
    Hzx(srcx,bound_width:ydim-bound_width)=Hzx(srcx,bound_width:ydim-bound_width)-dt./mu0./delta.*Einc;
    Hzy(srcx,bound_width:ydim-bound_width)=Hzy(srcx,bound_width:ydim-bound_width)-dt./mu0./delta.*Einc;
    
end

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEDIUM PROPS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% DO DIFFUSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_Diffusion()
    global T_e T_i Pinst 
    global T_e_old T_i_old
    global n1 n2 n11 n21
    global dt delta
    global material
        
    Delta_Tx=-2.*T_e(n1+1:n2,n11+1:n21) ... 
                + T_e(n1:n2-1, n11+1:n21)  + ~material(n1:n2-1,n11+1:n21)   .*T_e(n1+1:n2,n11+1:n21) ... %left
                + T_e(n1+2:n2+1,n11+1:n21) + ~material(n1+2:n2+1,n11+1:n21) .*T_e(n1+1:n2,n11+1:n21); %right

    Delta_Ty=-2.*T_e(n1+1:n2,n11+1:n21) ...
                 +T_e(n1+1:n2,n11:n21-1)   + ~material(n1+1:n2,n11:n21-1)   .* T_e(n1+1:n2,n11+1:n21) ... %down
                 +T_e(n1+1:n2,n11+2:n21+1) + ~material(n1+1:n2,n11+2:n21+1) .* T_e(n1+1:n2,n11+1:n21); %up
    
    foo = 1e-4.*dt./delta.^2; % temperature-conduction = kappa/Cp ~ 100 mm^2/s (most metals)
    gamma = 5.69e17; % electron-phonon coupling in W/K/m^3
    Cp = 40350; % volumetric heat capacity in J/K/m^3    
    %Cp is approx linear: Cp=134.5*T_e
    
    %electron-phonon energy-coupling accounted for by Delta_e_ph
    %heat diffusion of ions/lattice can be neglected on this time-scale
    %because the lattice heat conductivity is very low compared to electron
    %gas for metals
    %Note, that all the thermophysical as well as the optical
    %parameters are dependent on a variety of
    %properties such as Temperature, density, degree of ionization.
    %This model is very likely to produce unrealistic estimates beyond the
    %Fermi-Temperature (~100,000 K)
    %
    Delta_e_ph = gamma.*dt./Cp.*(T_e(n1+1:n2,n11+1:n21)-T_i(n1+1:n2,n11+1:n21));
    
    %Note: I do not check the cfl-criterion for heat diffusion here, 
    %because the diffusion step is in sync with the maxwell-solver and
    %the FDTD-timestep is very likely much lower than what
    %this criterion would require anyway
    %
    T_e(n1+1:n2,n11+1:n21) = T_e(n1+1:n2,n11+1:n21) + foo.*(Delta_Tx+Delta_Ty) ... 
                            +Pinst(n1+1:n2,n11+1:n21).*dt./Cp ...
                            -Delta_e_ph;
                        
    T_i(n1+1:n2,n11+1:n21) = T_i(n1+1:n2,n11+1:n21) + Delta_e_ph;
                        
    
end