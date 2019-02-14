function [ Bz,BR ] = BmagnmirnvMulti( Z_filament,R_filament,I_filament,r_mirnv,z_mirnv)

%-------------------------------------------------------------------------%
%---Filament is in the R-Y plane and Magnetic Field is Evaluated -------------%
%-------------at the Mirnv coordinate in the R-Z plane------------------------%
%-------------------------------------------------------------------------%

Zc=Z_filament;
I=I_filament;

turns=1; %%%% porque é modelo de filamentos

N=100;   % No of grids in the coil ( X-Y plane)
u0=4*pi*0.001;   % [microWb/(A cm)]
phi=0:2*pi/(N-1):2*pi; % For describing a circle (coil)

Rc=R_filament*cos(phi); % R-coordinates of the filament
Yc=R_filament*sin(phi); % Y-coordinates of the filament
% Zc(1:25)=Zc;

%Lets obtain the position vectors from dl to the mirnov 
%%% mirnov is localized in the plane (y=0)
for i=1:N-1
   RR(i)=r_mirnv-0.5*(Rc(i)+Rc(i+1)); 
   Rz(i)=z_mirnv-Zc;
   Ry(i)=-0.5*(Yc(i)+Yc(i+1));     
    dlR(i)=Rc(i+1)-Rc(i);
    dly(i)=Yc(i+1)-Yc(i);
end
RR(N)=r_mirnv-0.5*(Rc(N)+Rc(1)); 
   Rz(N)=z_mirnv-Zc;
   Ry(N)=-0.5*(Rc(N)+Rc(1));     
    dlR(N)=-Rc(N)+Rc(1);
    dly(N)=-Yc(N)+Yc(1);

%%dl x r
for i=1:N
Rcross(i)=-dly(i).*Rz(i);
Ycross(i)=dlR(i).*Rz(i);
Zcross(i)=(dly(i).*RR(i))-(dlR(i).*Ry(i));
R(i)=sqrt(RR(i).^2+Rz(i).^2+Ry(i).^2);
end


%dB=m0/4pi (Idl x r)/r^2 
BR1=(I*u0./(4*pi*(R.^3))).*Rcross;
Bz1=(I*u0./(4*pi*(R.^3))).*Zcross;
By1=(I*u0./(4*pi*(R.^3))).*Ycross;

BR=0;       % Initialize sum magnetic field to be zero first
Bz=0;
By=0;

for i=1:N   % loop over all current elements along coil    
    BR=BR+BR1(i);
    Bz=Bz+Bz1(i);
    By=By+By1(i);
end

BR=turns*BR;
Bz=turns*Bz;
By=turns*By; %%% units=[uWb / cm^2]

vector=[Z_filament-z_mirnv,R_filament-r_mirnv];%Vector from center of chamber to mirnov center
unit_vec=[vector]./norm(vector); %% Unit vector
norm_vec=[-unit_vec(1),unit_vec(2)];%%%  Normal vector, coil direction
 Bmirn=abs(BR*unit_vec(1)+Bz*unit_vec(2));
 Bmirn=sqrt(BR^2+Bz^2+By^2);
 Bmirn=dot([BR,Bz],norm_vec);
Bmirn=0.01*Bmirn;%fator de 0.01 pra ter [T] 
end