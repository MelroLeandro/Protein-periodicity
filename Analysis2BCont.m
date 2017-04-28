%%
% Analysing window 1
%

clear all
clc
close all
format long

addpath pdbTools CircStat

KC_C = 326;
KN_C = 250;

MeandC_CA= 1.52;
MeandN_CA= 1.45;

bin=pi/180*8.6;
kappa=5;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

%filter=@(dC_CA,dN_CA,dP_plane) not(abs(dC_CA-1.99)>0.5 || abs(dN_CA-1.78)>0.4 || abs(dP_plane-1.47)>0.3);
filter=@(dC_CA,dN_CA,dP_plane) not(abs(dC_CA-1.99)>1.18 || abs(dN_CA-1.78)>0.91 || abs(dP_plane-1.47)>0.55);
%filter=@(dC_CA,dN_CA,dP_plane) not(dC_CA<2*1.76 || dN_CA<2*1.7 || dP_plane<2*1.7);
%filter=@(dC_CA,dN_CA,dP_plane)1;

%aminos={'PHE'; 'ASP'; 'THR'; 'ARG'; 'TRP'; 'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
aminos1={'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
aminos2={'PHE'; 'ASP'; 'THR'; 'ARG'; 'TRP'; 'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
varName={'phi', 'psi', 'dC_CA','dN_CA','dN_C','dP_plane', 'ang_N_CA_C','R1','R2','R3','R4'};
for j=1:length(aminos1) 
    for i=1:length(aminos2)   
        aminoName1=aminos1{j}
        aminoName2=aminos2{i}
        Data=importdata(['aminoData/' aminoName1 aminoName2 '.cvs']);

        if ~isempty(Data)

            [l,col]=size(Data);

            data=removenan2(Data,filter,l); % clean and select data

            %%
            % Amino discrete characteristic
            %
            phi1= data(:,1);
            psi1= data(:,2);

            % bounds
            C_CA1= data(:,3);
            N_CA1= data(:,4);
            dplane1= data(:,6);
            
            phi2= data(:,8);
            psi2= data(:,9);
            
            omega= data(:,10);

            % bounds
            C_CA2= data(:,11);
            N_CA2= data(:,12);
            dplane2= data(:,14);
            
            

            [f,w1]=circ_vmpar(phi1);
            [f,w2]=circ_vmpar(psi1);
            [f,w3]=circ_vmpar(phi2);
            [f,w4]=circ_vmpar(psi2);
            [f,w5]=circ_vmpar(omega);
            %lambda = circ_density( phi,psi, w1,w2);
            
            %% Deep 1
            [p nbins]=circ_vmpdf2(180, phi1, phi2, w1, w3);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(1,:,:)=-log(p);

            title([aminoName1 aminoName2 ': phi1 vs phi2']);
            xlabel('phi1');
            ylabel('phi2');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_1',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, phi1, phi2, w1, w3);
            [d2 nbins] = circ_d2vmpdf2(180, phi1, phi2, w1, w3);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(1,:,:)=-log(d1p);
            aminoCont.dfy(1,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) phi1 vs phi2']);
            xlabel('phi1');
            ylabel('phi2');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad1',aminoName1,aminoName2,'.eps'))


            
            %% Deep 2
            [p nbins]=circ_vmpdf2(180, psi1, psi2, w2, w4);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(2,:,:)=-log(p);

            title([aminoName1 aminoName2 ': psi1 vs psi2']);
            xlabel('psi1');
            ylabel('psi2');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_2',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, psi1, psi2, w2, w4);
            [d2 nbins] = circ_d2vmpdf2(180, psi1, psi2, w2, w4);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(2,:,:)=-log(d1p);
            aminoCont.dfy(2,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) psi1 vs psi2']);
            xlabel('psi1');
            ylabel('psi2');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad1',aminoName1,aminoName2,'.eps'))

            
            
            %% Deep 3
            [p nbins]=circ_vmpdf2(180, psi1, phi2, w2, w3);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(3,:,:)=-log(p);

            title([aminoName1 aminoName2 ': psi1 vs phi2']);
            xlabel('psi1');
            ylabel('phi2');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_Grad2',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, psi1, phi2, w2, w3);
            [d2 nbins] = circ_d2vmpdf2(180, psi1, phi2, w2, w3);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(3,:,:)=-log(d1p);
            aminoCont.dfy(3,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) psi1 vs phi2']);
            xlabel('psi1');
            ylabel('phi2');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad2',aminoName1,aminoName2,'.eps'))

            
            
            %% Deep 4
            [p nbins]=circ_vmpdf2(180, phi1, psi2, w1, w4);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(4,:,:)=-log(p);

            title([aminoName1 aminoName2 ': phi1 vs psi2']);
            xlabel('phi1');
            ylabel('psi2');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_4',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, phi1, psi2, w1, w4);
            [d2 nbins] = circ_d2vmpdf2(180, phi1, psi2, w1, w4);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(4,:,:)=-log(d1p);
            aminoCont.dfy(4,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) phi1 vs psi2']);
            xlabel('phi1');
            ylabel('psi2');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad4',aminoName1,aminoName2,'.eps'))
            
            
            %%
            % Deep 5 opega 1
            
            [p nbins]=circ_vmpdf2(180, psi1, omega, w2, w5);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(4,:,:)=-log(p);

            title([aminoName1 aminoName2 ': psi1 vs omega']);
            xlabel('psi1');
            ylabel('omega');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_5',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, psi1, omega, w2, w5);
            [d2 nbins] = circ_d2vmpdf2(180, psi1, omega, w2, w5);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(4,:,:)=-log(d1p);
            aminoCont.dfy(4,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) psi1 vs omega']);
            xlabel('psi1');
            ylabel('omega');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad5',aminoName1,aminoName2,'.eps'))
            
            %%
            %% Deep 6 opega 2
            
            [p nbins]=circ_vmpdf2(180, phi2, omega, w3, w5);
            
            
            %figure
            imagesc(-log(p));

            %pause

            aminoCont.f(4,:,:)=-log(p);

            title([aminoName1 aminoName2 ': phi2 vs omega']);
            xlabel('phi2');
            ylabel('omega');
            %%axis([-180 180 -180 180]);
            colorbar;

            print(gcf,'-dpsc2',strcat('images/Rama_6',aminoName1,aminoName2,'.eps'))


            [d1 nbins] = circ_d1vmpdf2(180, phi2, omega, w3, w5);
            [d2 nbins] = circ_d2vmpdf2(180, phi2, omega, w3, w5);
            d1p=-d1./p;
            d2p=-d2./p;
            aminoCont.dfx(4,:,:)=-log(d1p);
            aminoCont.dfy(4,:,:)=-log(d2p);

            %figure
            imagesc(sqrt(d1p.^2+d2p.^2));
            title([aminoName1 aminoName2 ':norm(Grad) phi2 vs omega']);
            xlabel('phi2');
            ylabel('omega');
            %%axis([-180 180 -180 180]);
            colorbar;


            print(gcf,'-dpsc2',strcat('images/Rama_Grad6',aminoName1,aminoName2,'.eps'))
            
            %%
            
            %

            save(['contData/' aminoName1 aminoName2 '.mat'],'aminoCont');

        end


    end
end


