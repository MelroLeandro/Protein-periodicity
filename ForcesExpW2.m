function aminoContW2 = ForcesExpW2( aminoNames, binAng)
%
% Computes dihedral and rotamer angles for a given amino 
%

addpath pdbTools CircStat

SbinAng=int2str(binAng);
bin=pi/180*binAng;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(abs(dC_CA-1.52)>0.03 || abs(dN_CA-1.47)>0.03 || abs(dP_plane-1.32)>0.03 || Bfactor>40);

Data=importdata(['aminoData/' aminoNames '.cvs']);

if ~isempty(Data)

[l,col]=size(Data);

%data=removenan(Data,filter); % clean and select data

%%
% Amino discrete characteristic
%
phi1= Data(:,1);
psi1= Data(:,2);


phi2= Data(:,8);
psi2= Data(:,9);

Temp= Data(:,10);


fprintf('\t\t\t\t Prossecing phi1,psi2\n');
[p nbins total]=circ_exp(binAng,  phi1, psi2);

aminoContW2.f(1,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
[d1 d2 nbins] = circ_d1exp(binAng, phi1,psi2,total);

d1p=-d1./p;
d2p=-d2./p;
aminoContW2.dfx(1,:,:)=d1p;
aminoContW2.dfy(1,:,:)=d2p;

%%

fprintf('\t\t\t\t Prossecing psi1,phi2\n');
[p nbins total]=circ_exp(binAng,  psi1, phi2);

aminoContW2.f(2,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
[d1 d2 nbins] = circ_d1exp(binAng,psi1,phi2,total);
d1p=-d1./p;
d2p=-d2./p;
aminoContW2.dfx(2,:,:)=d1p;
aminoContW2.dfy(2,:,:)=d2p;

%%

fprintf('\t\t\t\t Prossecing psi1,phi1\n');
[p nbins total]=circ_exp(binAng,  psi1, phi1);

aminoContW2.f(3,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi1)\n');
[d1 d2 nbins] = circ_d1exp(binAng,psi1,phi1,total);

d1p=-d1./p;
d2p=-d2./p;
aminoContW2.dfx(3,:,:)=d1p;
aminoContW2.dfy(3,:,:)=d2p;

save(['contData/Exp' aminoNames SbinAng '.mat'],'aminoContW2');        
end


