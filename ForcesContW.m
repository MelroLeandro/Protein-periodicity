function aminoContW = ForcesContW2( aminoNames, binAng)
%
% Computes dihedral and rotamer angles for a given amino 
%

addpath pdbTools CircStat


bin=pi/180*binAng;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(abs(dC_CA-1.52)>0.03 || abs(dN_CA-1.47)>0.03 || abs(dP_plane-1.32)>0.03 || Bfactor>40);

Data=importdata(['aminoDataOld/' aminoNames '.cvs']);

if ~isempty(Data)

[l,col]=size(Data);

%data=removenan(Data,filter); % clean and select data

%%
% Amino discrete characteristic
%
phi1= Data(:,1);
psi1= Data(:,2);

% bounds
C_CA1= Data(:,3);
N_CA1= Data(:,4);
dplane1= Data(:,5);

phi2= Data(:,8);
psi2= Data(:,9);

Temp= Data(:,10);

% bounds
C_CA2= Data(:,11);
N_CA2= Data(:,12);
dplane2= Data(:,13);

[f,w1]=circ_vmpar(phi1);
[f,w2]=circ_vmpar(psi2);

fprintf('\t\t\t\t Prossecing phi1,psi2\n');
[p nbins]=circ_vmpdf2(binAng, phi1, psi2, w1, w2);
[p nbins]=circ_vmpdf2(binAng, phi1, psi2, w1, w2,p);

aminoContW2.f(1,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, phi1,psi2, w1,w2);
%[d1 nbins] = circ_d1vmpdf2(binAng, phi1,psi2, w1,w2,p);

%[d2 nbins] = circ_d2vmpdf2(binAng, phi1,psi2, w1,w2);
%[d2 nbins] = circ_d2vmpdf2(binAng, phi1,psi2, w1,w2,p);

[d1p,d2p] = imgradientxy(-log(p));
aminoContW2.dfx(1,:,:)=d1p;
aminoContW2.dfy(1,:,:)=d2p;

[f,w1]=circ_vmpar(phi1);
[f,w2]=circ_vmpar(psi2);

%%

fprintf('\t\t\t\t Prossecing psi1,phi2\n');
[p nbins]=circ_vmpdf2(binAng, psi1, phi2, w1, w2);
[p nbins]=circ_vmpdf2(binAng, psi1, phi2, w1, w2,p);

aminoContW2.f(2,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, psi1,phi2, w1,w2);
%[d1 nbins] = circ_d1vmpdf2(binAng, psi1,phi2, w1,w2,p);

%[d2 nbins] = circ_d2vmpdf2(binAng, psi1,phi2, w1,w2);
%[d2 nbins] = circ_d2vmpdf2(binAng, psi1,phi2, w1,w2,p);

[d1p,d2p] = imgradientxy(-log(p));
aminoContW2.dfx(2,:,:)=d1p;
aminoContW2.dfy(2,:,:)=d2p;

%%

fprintf('\t\t\t\t Prossecing psi1,phi1\n');
[p nbins]=circ_vmpdf2(binAng, psi1, phi1, w1, w2);
[p nbins]=circ_vmpdf2(binAng, psi1, phi1, w1, w2,p);

aminoContW2.f(3,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi1)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, psi1,phi1, w1,w2);
%[d1 nbins] = circ_d1vmpdf2(binAng, psi1,phi1, w1,w2,p);

%[d2 nbins] = circ_d2vmpdf2(binAng, psi1,phi1, w1,w2);
%[d2 nbins] = circ_d2vmpdf2(binAng, psi1,phi1, w1,w2,p);

[d1p,d2p] = imgradientxy(-log(p));
aminoContW2.dfx(3,:,:)=d1p;
aminoContW2.dfy(3,:,:)=d2p;

save(['contData/Mises' aminoNames '.mat'],'aminoContW2');        
end


