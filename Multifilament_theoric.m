%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Plasma current centroid position reconstruction%%%%%%%%
%%%%%% Multifilaments,7 filaments , Pitonti's model %%%%%%%%%%
tic
close all
clear all
%%% Load Shot
load('shot_45520.mat');
time=1e-3*data.time; %%%% time in ms


%%% Draw the vessel
th = 0:pi/50:2*pi;
xvess = 9 * cos(th)+46;
yvess = 9 * sin(th) ;

%%% Mirnov positions
ang=-15;
for i=1:12
R_mirn(i)=9.35*cosd(ang)+46;
z_mirn(i)=9.35*sind(ang);
ang=ang-30;
end

%%%%%% Lets draw the plasma filaments
th1 = 0:pi/50:2*pi;
R_pls=46;
z_plsm= 0;

R_filaments(1)=46;
z_filaments(1)=0;
degr=0;
radius=4; %%% in [cm] (distance from the center of the chamber to the filaments)

for i=2:7
    R_filaments(i)=(46)+radius*cosd(degr);
    z_filaments(i)=radius*sind(degr);
    degr=degr+60;
end


%%%Experimental mesurements[Wb]

Mirnv_10_fact=1.2803;
time_ins=116;
time_index=find(time == 116); %%% Select a time moment where there is plasma current! in [ms]

%%%%%%%%%% Find the exprimental values for that time moment

%%%%without external flux correction

Mirnv_flux(:)=data.mirnv_corr(:,time_index);
Mirnv_flux(10)=Mirnv_10_fact*Mirnv_flux(10);

Mirnv_flux_corr(:)=data.mirnv_corr_flux(:,time_index);
Mirnv_flux_corr(10)=Mirnv_10_fact*Mirnv_flux_corr(10);

%%%%% Let's go from [Wb] to {T]
%%%% fmincon needs data to be double
Mirnv_B_exp=double(Mirnv_flux/(50*49e-6)); %%%% [T]
Mirnv_B_exp_corr=double(Mirnv_flux_corr/(50*49e-6)); %%%% [T]


 
%%%% Matrix whose elements gives the contribution  to the measuremnt i  to
%%%% a unitary current in the filament j [T]
for i=1:12
    for j=1:7
   
         Mfp(i,j)=Bmagnmirnv(z_filaments(j),R_filaments(j),1,R_mirn(i),z_mirn(i)) ;
    end
end

% Mfp=Mfp/mean(mean(Mfp));
Mpf=pinv(Mfp);
Mpf=Mpf/mean(mean(Mpf));

L=1e-6*[15.922
15.922
15.922
16.42
16.716
16.447
16.476
15.702
16.435
8.555
16.252
14.787];

I_filament=Mpf*(Mirnv_flux_corr'./L);
toc
