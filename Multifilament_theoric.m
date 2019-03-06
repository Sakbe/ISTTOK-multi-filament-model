%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Plasma current centroid position reconstruction%%%%%%%%
%%%%%% Multifilaments,7 filaments , SVD matrix %%%%%%%%%%
tic
close all
clear all
%%% Load Shot
load('shot_45860.mat');
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
degr=10;
radius=5.0; %%% in [cm] (distance from the center of the chamber to the filaments)
nfil=12; %%% Number of filaments
deg_fact=360/(nfil-1);

for i=2:nfil
    R_filaments(i)=(46)+radius*cosd(degr);
    z_filaments(i)=radius*sind(degr);
    degr=degr+deg_fact;
end


%%%Experimental mesurements[Wb]

%Mirnv_10_fact=1.2803;
time_ins=120;
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

%%%%% Minimization 7 degrees of freedom
%   fval_multi_corr=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,z_filaments(1),R_filaments(1),x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),R_filaments,z_filaments,R_mirn,z_mirn,nfil),[500,500,500,500,500,500,500,500,500,500,500])
%   fval_multi_corr=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,z_filaments(1),R_filaments(1),x,R_filaments,z_filaments,R_mirn,z_mirn,nfil),[500,500,500,500,500,500,500,500,500,500,500])


 %%
% xx_multi_corr=BmagnMultiModule_correct(z_filaments(1),R_filaments(1),fval_multi_corr,R_filaments,z_filaments,R_mirn,z_mirn,nfil);
%%%% Matrix whose elements gives the contribution  to the measuremnt i  to
%%%% a unitary current in the filament j [T]
for i=1:12
    for j=1:nfil
   
         Mfp(i,j)=Bmagnmirnv(z_filaments(j),R_filaments(j),1,R_mirn(i),z_mirn(i)) ;
    end
end

Mpf=pinv(Mfp);
I_filament=Mpf*(Mirnv_B_exp_corr');

%% Calculate Biot-Savart with the current values from SVD decomposition
xx_multi_SVD=BmagnMultiModule_correct(z_filaments(1),R_filaments(1),I_filament,R_filaments,z_filaments,R_mirn,z_mirn,nfil);

% for i=1:12
%    xx_multi_theo(i) =0;
%     for j=1:7
%         
% xx_multi_theo(i)=Bmagnmirnv(z_filaments(j),R_filaments(j),I_filament(j),R_mirn(i),z_mirn(i)) +xx_multi_theo(i);
%     end
% end
%% Error
RMSE_optim_theo=sqrt(mean((xx_multi_SVD(:)-Mirnv_B_exp_corr(:)).^2));



%% Compute Centroid position
I_filament_all=Mpf*(data.mirnv_corr_flux)/(49*50*1e-6);
for(i=1:length(I_filament_all))
z0(i)=0.01*sum((z_filaments'.*I_filament_all(:,i)))./sum(I_filament_all(:,i));
r0(i)=0.01*sqrt(sum((R_filaments'.^2).*I_filament_all(:,i))./sum(I_filament_all(:,i)))-0.46;
sumIfil(i)=sum(I_filament_all(:,i));
end

for(i=1:length(z0))
if(imag(z0(i)~=0))
z0(i)=0.085;
end
if(z0(i)>0.085||z0(i) <-0.085)
z0(i)=0.085;
end
end

for(i=1:length(r0))
if(imag(r0(i))~=0)
r0(i)=0.085;
end
if(r0(i)>0.085||r0(i) <-0.085)
r0(i)=0.085;
end
end


%% Plotting
close all
figure(9)
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mirnv_B_exp_corr ,'-o')
hold on
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_SVD,'-s')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_corr,'-*')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mfp*I_filament,'-s')
grid on
title(['Shot #45410  t= ',num2str(time_ins), '  Ip= (Multifilament flux corrected)'])
legend('Experimental Data corrected','Multifilament SVD-Mpf','Fminsearch')
xlabel('Mirnov #')
ylabel('Optimization [mT]')
axis equal

figure(10)
plot(time,sumIfil)
hold on
plot(time,data.Ip_magn)
grid on
xlabel('Time [ms]')
ylabel('Current [A]')
legend('sum (I fil)','Plasma current')

%%
figure(3)
plot(xvess,yvess,'k','linewidth',2)
hold on
plot(46,0,'.m','MarkerSize',790)
plot(R_mirn,z_mirn,'sk','MarkerSize',17)

% plot(fval_multi_corr(2),fval_multi_corr(1),'.k','MarkerSize',20)
for i=1:nfil
    plot(R_filaments(i),z_filaments(i),'.b','MarkerSize',20)
end
    for i = 1:12
    text(R_mirn(i),z_mirn(i),num2str(i),'Color','r','FontSize',13) 


end
plot([35,57],[0,0])
plot([46,46],[-9,9],'k')
text(57,0,'LFS','FontSize',15)
text(33,0,'HFS','FontSize',15)
ylim([-11,11])
xlabel('R[cm]')
ylabel('Z[cm]')
grid on
axis equal
toc
return


figure(1)
plot(time,R0)
grid on
xlabel('Time [ms]')
ylabel('Position [m]')
title('Radial centroid position')

figure(2)
plot(time,z0)
grid on
xlabel('Time [ms]')
ylabel('Position [m]')
title('Vertical centroid position')