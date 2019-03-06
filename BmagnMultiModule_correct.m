function Bmirn=BmagnMultiModule_correct(Z_filament,R_filament,I_filaments,R_filaments,z_filaments,r_mirnv,z_mirnv,nfil)


for i=1:12
    for j=1:nfil
vector(i,j,[1,2])=[z_filaments(j)-z_mirnv(i),R_filaments(j)-r_mirnv(i)];%Vector from center of chamber to mirnov center
unit_vector(i,j,:)=[vector(i,j,1),vector(i,j,2)]./norm([vector(i,j,1),vector(i,j,2)]); %% Unit vector
norm_vec(i,j,:)=[unit_vector(i,j,2),-unit_vector(i,j,1)];%%%  Normal vector, coil direction
    end
end

for i=1:12  
[Bz(i,1),BR(i,1)]=BmagnmirnvMulti(Z_filament,R_filament,I_filaments(1),r_mirnv(i),z_mirnv(i));
for j=2:nfil
[Bz(i,j),BR(i,j)]=BmagnmirnvMulti(z_filaments(j),R_filaments(j),I_filaments(j),r_mirnv(i),z_mirnv(i));
end 
end




%%% Lets make the sum as it has to be done
for i=1:12
Bmirn(i)=0;
for j=1:nfil
    Bmirn(i)=Bmirn(i)+dot([Bz(i,j),BR(i,j)],[norm_vec(i,j,1),norm_vec(i,j,2)]);
end
end

Bmirn=0.01*Bmirn;%fator de 0.01 pra ter [T] 



end