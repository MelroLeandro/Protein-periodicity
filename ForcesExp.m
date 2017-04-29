function aminoCont = ForcesExp( aminoName, binAng)
%
% Computes dihedral and rotamer angles for a given amino 
%

addpath pdbTools CircStat

SbinAng=int2str(binAng);
bin=pi/180*binAng;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(abs(dC_CA-1.52)>0.03 || abs(dN_CA-1.47)>0.03 || abs(dP_plane-1.32)>0.03 || Bfactor>40);

    
Data=importdata(['aminoData/' aminoName '.cvs']);

if ~isempty(Data)

[l,col]=size(Data);

data=removenan(Data,filter); % clean and select data

%%
% Amino discrete characteristic
%
phi= data(:,1);
psi= data(:,2);


fprintf('\t\t\t\t Prossecing phi,psi\n');
[p nbins total]=circ_exp(binAng, phi, psi);

aminoCont.f(1,:,:)=-log(p);

fprintf('\t\t\t\t Prossecing Grad(phi,psi)\n');
[d1 d2 nbins] = circ_d1exp(binAng, phi,psi,total);

d1p=-d1./p;
d2p=-d2./p;
aminoCont.dfx(1,:,:)=d1p;
aminoCont.dfy(1,:,:)=d2p;

aminoCont.N_rota=col-8;

if col> 8
            chsi1= data(:,9);
            
            fprintf('\t\t\t\t Prossecing phi,chsi1\n');
            
            [p nbins total]=circ_exp(binAng,  phi, chsi1);
            aminoCont.f(2,:,:)=-log(p);
            
            fprintf('\t\t\t\t Prossecing Grad(phi,chsi1)\n');
            [d1 d2 nbins] = circ_d1exp(binAng, phi, chsi1,total);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(2,:,:)=d1p;
            aminoCont.dfy(2,:,:)=d2p;
            
            fprintf('\t\t\t\t Prossecing psi,chsi1\n');
            [p nbins total]=circ_exp(binAng,  psi, chsi1);
            aminoCont.f(3,:,:)=-log(p);
            
            fprintf('\t\t\t\t Prossecing Grad(psi,chsi1)\n');
            [d1 d2 nbins] = circ_d1exp(binAng, psi, chsi1,total);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(3,:,:)=d1p;
            aminoCont.dfy(3,:,:)=d2p;


            if col > 9                        
                chsi2= data(:,10);

                fprintf('\t\t\t\t Prossecing chsi1,chsi2\n');
                [p nbins total]=circ_exp(binAng, chsi1, chsi2);
                aminoCont.f(4,:,:)=-log(p);
                
                fprintf('\t\t\t\t Prossecing Grad(chsi1,chsi2)\n');
                
                [d1 d2 nbins] = circ_d1exp(binAng, chsi1, chsi2,total);
                d1p=-d1./p;
                d2p=-d2./p;
                aminoCont.dfx(4,:,:)=d1p;
                aminoCont.dfy(4,:,:)=d2p;

                if col > 10

                    chsi3= data(:,11);
                    
                    fprintf('\t\t\t\t Prossecing chsi2,chsi3\n');
                    [p nbins total]=circ_exp(binAng, chsi2, chsi3);
                    aminoCont.f(5,:,:)=-log(p);

                    fprintf('\t\t\t\t Prossecing Grad(chsi2,chsi3)\n');
                    [d1 d2 nbins] = circ_d1exp(binAng, chsi2, chsi3,total);
                    d1p=-d1./p;
                    d2p=-d2./p;
                    aminoCont.dfx(5,:,:)=d1p;
                    aminoCont.dfy(5,:,:)=d2p;

                    if col > 11

                        chsi4= data(:,12);
                        
                        fprintf('\t\t\t\t Prossecing chsi4,chsi5\n');
                        [p nbins total]=circ_exp(binAng, chsi3, chsi4);
                        aminoCont.f(6,:,:)=-log(p);
                        fprintf('\t\t\t\t Prossecing Grad(chsi4,chsi5)\n');
                        [d1 d2 nbins] = circ_d1exp(binAng, chsi3, chsi4,total);
                        d1p=-d1./p;
                        d2p=-d2./p;
                        aminoCont.dfx(6,:,:)=d1p;
                        aminoCont.dfy(6,:,:)=d2p;

                    end
                end
            end
end


save(['contData/Exp' aminoName SbinAng '.mat'],'aminoCont');        
end


