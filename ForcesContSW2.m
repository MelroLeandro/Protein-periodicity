function aminoContW2 = ForcesContSW2( aminoNames, binAng)
%
% Computes dihedral and rotamer angles for a given amino 
%

addpath pdbTools CircStat

SbinAng=int2str(binAng);

bin=pi/180*binAng;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

Data=importdata(['aminoData/' aminoNames '.cvs']);

if ~isempty(Data)

[l,col]=size(Data);

data=removenan2(Data); % clean and select data

%%
% Amino discrete characteristic
%
phi1= data(:,1);
psi1= data(:,2);

% bounds
C_CA1= data(:,3);
N_CA1= data(:,4);
dplane1= data(:,5);

phi2= data(:,8);
psi2= data(:,9);

Temp= data(:,10);

% bounds
C_CA2= data(:,11);
N_CA2= data(:,12);
dplane2= data(:,13);

[f,wphi1]=circ_vmpar(phi1);
[f,wpsi1]=circ_vmpar(psi1);
[f,wphi2]=circ_vmpar(phi2);
[f,wpsi2]=circ_vmpar(psi2);

fprintf('\t\t\t\t Prossecing phi1,psi2\n');
[p nbins]=circ_vmpdf2(binAng, phi1, psi2, wphi1, wpsi2);

aminoContW2.f(1,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, phi1,psi2, wphi1,wpsi2);

%[d2 nbins] = circ_d2vmpdf2(binAng, phi1,psi2, wphi1,wpsi2);

%d1p=-d1./p;
%d2p=-d2./p;
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));

aminoContW2.dfx(1,:,:)=gather(d1p);
aminoContW2.dfy(1,:,:)=gather(d2p);

%%

fprintf('\t\t\t\t Prossecing phi1,phi2\n');
[p nbins]=circ_vmpdf2(binAng, phi1, phi2, wpsi1, wphi2);

aminoContW2.f(2,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,phi2)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, phi1,phi2, wpsi1,wphi2);
%
%[d2 nbins] = circ_d2vmpdf2(binAng, phi1,phi2, wpsi1,wphi2);

%d1p=-d1./p;
%d2p=-d2./p;
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));

aminoContW2.dfx(2,:,:)=gather(d1p);
aminoContW2.dfy(2,:,:)=gather(d2p);

%%

fprintf('\t\t\t\t Prossecing psi1,phi2\n');
[p nbins]=circ_vmpdf2(binAng, psi1, phi2, wpsi1, wphi2);

aminoContW2.f(3,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(psi1,psi2)\n');
%[d1 nbins] = circ_d1vmpdf2(binAng, psi1,phi2, wpsi1,wphi2);
%
%[d2 nbins] = circ_d2vmpdf2(binAng, psi1,phi2, wpsi1,wphi2);

%d1p=-d1./p;
%d2p=-d2./p;
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));

aminoContW2.dfx(3,:,:)=gather(d1p);
aminoContW2.dfy(3,:,:)=gather(d2p);

save(['contData/MisesS' aminoNames SbinAng '.mat'],'aminoContW2');        
end


