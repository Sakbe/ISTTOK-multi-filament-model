function error=ErrorMirnFuncMultiFilam(Mirnv_B_exp,Z_filament,R_filament,I_filament,...
    I_filament2,I_filament3,I_filament4,I_filament5,I_filament6,I_filament7,R_filaments,z_filaments,r_mirnv,z_mirnv)


for i=1:12
vector(i,[1,2])=[z_filaments(1)-z_mirnv(i),R_filaments(1)-r_mirnv(i)];%Vector from center of chamber to mirnov center
unit_vec(i,[1,2])=[vector(i,[1,2])]./norm(vector(i,[1,2])); %% Unit vector
norm_vec(i,[1,2])=[unit_vec(i,2),-unit_vec(i,1)];%%%  Normal vector, coil direction
end

%%%% Create vector of the current from the filaments
I_filaments=[I_filament,I_filament2,I_filament3,I_filament4,I_filament5,I_filament6,I_filament7];

%%%%% Calculate the contribution of each filament to each mirnov coil
for i=1:12    
[Bz(i,1),BR(i,1)]=BmagnmirnvMulti(Z_filament,R_filament,I_filaments(1),r_mirnv(i),z_mirnv(i));
for j=2:7
[Bz(i,j),BR(i,j)]=BmagnmirnvMulti(z_filaments(j),R_filaments(j),I_filaments(j),r_mirnv(i),z_mirnv(i));
end  
end

%%%%% Vectorial sum from the effect of the filaments on each mirnov coil
for i=1:12
Bz_tot(i)=0;
BR_tot(i)=0;
end

for i=1:12
for j=1:7
     BR_tot(i)=BR_tot(i)+BR(i,j);
     Bz_tot(i)=Bz_tot(i)+Bz(i,j);
end
end

%%% Calculate the projection 
for i=1:12
Bmirn(i)=dot([Bz_tot(i),BR_tot(i)],norm_vec(i,[1,2]));
end

Bmirn=0.01*Bmirn;%fator de 0.01 pra ter [T] 


error=sum(abs(Mirnv_B_exp-Bmirn));

end