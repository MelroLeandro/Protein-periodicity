function aminoContW2 = ForcesContSW3( aminoName1, aminoName2, aminoName3, binAng)
%
% Computes dihedral and rotamer angles for a given amino 
%

addpath pdbTools CircStat


bin=pi/180*binAng;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

%filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(abs(dC_CA-1.52)>0.03 || abs(dN_CA-1.47)>0.03 || abs(dP_plane-1.32)>0.03 || Bfactor>40);
filter=@(dC_CA,dN_CA,dP_plane,Bfactor) 1;

%Data=fopen('/home/cortex/Desktop/Protein_measure/rotamer_w3.cvs','r');
file=fopen('C:\Users\Administrator\Desktop\protein_measure2\rotamer_w3.cvs','r');

tline = fgetl(file);

Data=[];
while ischar(tline)
    %fprintf('%s\n',tline(1:11));
    if strcmp(tline(1:11),[aminoName1 ',' aminoName2 ',' aminoName3])
        Data=[Data tline];
    end
    tline = fgetl(file);
end

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

aminoContW2.f(1,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));
aminoContW2.dfx(1,:,:)=gather(d1p);
aminoContW2.dfy(1,:,:)=gather(d2p);

[f,w1]=circ_vmpar(phi1);
[f,w2]=circ_vmpar(psi2);

%%

fprintf('\t\t\t\t Prossecing psi1,phi2\n');
[p nbins]=circ_vmpdf2(binAng, psi1, phi2, w1, w2);

aminoContW2.f(2,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi2)\n');
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));
aminoContW2.dfx(2,:,:)=gather(d1p);
aminoContW2.dfy(2,:,:)=gather(d2p);

%%

fprintf('\t\t\t\t Prossecing psi1,phi1\n');
[p nbins]=circ_vmpdf2(binAng, psi1, phi1, w1, w2);

aminoContW2.f(3,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi1,psi1)\n');
I = gpuArray(p);
[d1p,d2p] = imgradientxy(-log(I));
aminoContW2.dfx(3,:,:)=gather(d1p);
aminoContW2.dfy(3,:,:)=gather(d2p);

save(['contData/MisesS',  aminoName1, aminoName2, aminoName3, '.mat'],'aminoContW2');        
end


