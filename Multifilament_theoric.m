%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Plasma current centroid position reconstruction%%%%%%%%
%%%%%% Multifilaments,7 filaments , SVD matrix %%%%%%%%%%
tic
close all
clear all
%%% Load Shot
load('shot_45825.mat');
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
radius=5; %%% in [cm] (distance from the center of the chamber to the filaments)

for i=2:7
    R_filaments(i)=(46)+radius*cosd(degr);
    z_filaments(i)=radius*sind(degr);
    degr=degr+60;
end


%%%Experimental mesurements[Wb]

%Mirnv_10_fact=1.2803;
time_ins=115;
time_index=find(time == time_ins); %%% Select a time moment where there is plasma current! in [ms]

%%%%%%%%%% Find the exprimental values for that time moment

%%%%without external flux correction

Mirnv_flux(:)=data.mirnv_corr(:,time_index);
%Mirnv_flux(10)=Mirnv_10_fact*Mirnv_flux(10);

Mirnv_flux_corr(:)=data.mirnv_corr_flux(:,time_index);
%Mirnv_flux_corr(10)=Mirnv_10_fact*Mirnv_flux_corr(10);

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

Mpf=pinv(Mfp);
I_filament=Mpf*(Mirnv_B_exp_corr');

%% Calculate Biot-Savart with the current values from SVD decomposition
xx_multi_SVD=BmagnMultiModule_correct(z_filaments(1),R_filaments(1),I_filament,R_filaments,z_filaments,R_mirn,z_mirn);

% for i=1:12
%    xx_multi_theo(i) =0;
%     for j=1:7
%         
% xx_multi_theo(i)=Bmagnmirnv(z_filaments(j),R_filaments(j),I_filament(j),R_mirn(i),z_mirn(i)) +xx_multi_theo(i);
%     end
% end
%% Error

RMSE_optim_theo=sqrt(mean((xx_multi_SVD(:)-Mirnv_B_exp_corr(:)))^2);


close all
figure(9)
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mirnv_B_exp_corr ,'-o')
hold on
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_SVD,'-s')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mfp*I_filament,'-s')
grid on
title(['Shot #45410  t= ',num2str(time_ins), '  Ip= (Multifilament flux corrected)'])
legend('Experimental Data corrected','Multifilament SVD-Mpf')
xlabel('Mirnov #')
ylabel('Optimization [mT]')
axis equal

toc
