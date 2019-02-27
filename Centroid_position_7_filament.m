%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Plasma current centroid position reconstruction%%%%%%%%
%%%%%% Multifilaments,7 filaments, 9 freedom degrees%%%%%%%%%%
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
radius=3.5; %%% in [cm] (distance from the center of the chamber to the filaments)

for i=2:7
    R_filaments(i)=(46)+radius*cosd(degr);
    z_filaments(i)=radius*sind(degr);
    degr=degr+60;
end


%%%Experimental mesurements[Wb]

Mirnv_10_fact=1.2803;
time_ins=111;
time_index=find(time == time_ins); %%% Select a time moment where there is plasma current! in [ms]

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



%%%% Optimization function, 7 filaments, 9 degrees of freedom
%%%%% Central filament - 3 dregrees of freedom (z,R,I)
%%%%%% 6 sorrounding filaments - 1 degree of freedom (I)

%%%%Lets put boundaries
if Mirnv_B_exp_corr(1)>0
low_bnd=[-2,43,-4000,-4000,-4000,-4000,-4000,-4000,-4000];
high_bnd=[2,50,0,0,0,0,0,0,0];
Ini_cond=[0.5,46.5,-1000,-500,-500,-500,-500,-500,-500];
else
    low_bnd=[-2,43,0,0,0,0,0,0,0];
high_bnd=[2,50,4000,4000,4000,4000,4000,4000,4000];
Ini_cond=[0.5,46.5,1000,500,500,500,500,500,500];
end

A=-eye(9);
A(1:2,1:9)=0;
b=zeros(9,1);

% fval_multi=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),[0.5,46.5,500,500,500,500,500,500,500])

fval_multi=fmincon(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),Ini_cond,[],[],[],[],low_bnd,high_bnd)

% fval_multi=fmincon(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),[0.5,46.5,1000,500,500,500,500,500,500],A,b)

%%%Externa fluxes corrected

% fval_multi_corr=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),[0.5,46.5,500,500,500,500,500,500,500])

fval_multi_corr=fmincon(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),Ini_cond,[],[],[],[],low_bnd,high_bnd)

% fval_multi_corr=fmincon(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),R_filaments,z_filaments,R_mirn,z_mirn),[0.5,46.5,1000,500,500,500,500,500,500],A,b)
 
%%%%Lets check how close is our minimization values to the experimental
%%%%ones by applaying Biot-Savart with them 

xx_multi=BmagnMultiModule_correct(fval_multi(1),fval_multi(2),fval_multi(3:9),R_filaments,z_filaments,R_mirn,z_mirn);

xx_multi_corr=BmagnMultiModule_correct(fval_multi_corr(1),fval_multi_corr(2),fval_multi_corr(3:9),R_filaments,z_filaments,R_mirn,z_mirn);


%%%% Error

RMSE_multi=sqrt(mean((xx_multi(:)-Mirnv_B_exp(:)))^2);
RMSE_multi_corr=sqrt(mean((xx_multi_corr(:)-Mirnv_B_exp_corr(:)))^2);
Error_multi=sum(abs(xx_multi-Mirnv_B_exp))/12;
Error_multi_corr=sum(abs(xx_multi_corr-Mirnv_B_exp_corr))/12;

%%%%%%%%%%Plotting
%%%%%%Multifilament plots


figure(8)
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mirnv_B_exp,'-o')
hold on
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi,'-*')
grid on
title('Shot #45410  t=195[ms]  Ip~4.1[kA] (Multifilament)')
legend('Experimental Data ','Biot-savart  (optimized )')
xlabel('Mirnov #')
ylabel('Optimization [mT]')
axis equal

figure(9)
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mirnv_B_exp_corr,'-o')
hold on
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_corr,'-*')
grid on
title(['Shot #45410  t= ',num2str(time_ins), '  Ip~4.1[kA] (Multifilament flux corrected)'])
legend('Experimental Data corrected','Biot-savart  (optimized )')
xlabel('Mirnov #')
ylabel('Optimization [mT]')
axis equal


%% 

%%%%%% Plasma, vessel and mirnov coil plot



figure(3)
plot(xvess,yvess,'k','linewidth',2)
hold on
plot(46,0,'.m','MarkerSize',790)
plot(R_mirn,z_mirn,'sk','MarkerSize',17)

plot(fval_multi_corr(2),fval_multi_corr(1),'.k','MarkerSize',20)
for i=2:7
    plot(R_filaments(i),z_filaments(i),'.b','MarkerSize',20)
end
    for i = 1:12
    text(R_mirn(i),z_mirn(i),num2str(i),'Color','r','FontSize',13) 


end

text(57,0,'LFS','FontSize',15)
text(33,0,'HFS','FontSize',15)
ylim([-11,11])
xlabel('R[cm]')
ylabel('Z[cm]')
grid on
axis equal
toc