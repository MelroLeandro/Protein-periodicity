%%
% Block of forces in each amino acid
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

aminos={'PHE'; 'ASP'; 'THR'; 'ARG'; 'TRP'; 'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
%aminos={'ARG'}
varName={'phi', 'psi', 'dC_CA','dN_CA','dN_C','dP_plane', 'ang_N_CA_C','R1','R2','R3','R4'};

for i=1:length(aminos)   
    aminoName=aminos{i}
    
    Data=importdata(['aminoData/' aminoName '.cvs']);
    
    if ~isempty(Data)
        
        [l,col]=size(Data);
        
        data=removenan(Data,filter,l); % clean and select data
        
        %%
        % Amino discrete characteristic
        %
        phi= data(:,1);
        psi= data(:,2);
        
        % bounds
        C_CA= data(:,3);
        N_CA= data(:,4);
        dplane= data(:,6);
        
        w=[];
        for ai=1:col
            [f,w1]=circ_vmpar(phi);
            w=[w w1]
        end
        
        M=circ_vmpdfBloc(180,data,w);
        %figure
        imagesc(-log(p));
        
        %pause
        
        aminoCont.f(1,:,:)=-log(p);
        
        title([aminoName ': phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
       
        print(gcf,'-dpsc2',strcat('images/Rama_',aminoName,'.eps'))
                    
        
        D = circ_d1vmpdf2(180,data,w);
        D=-D./M;
        
        aminoCont.dfx(1,:,:)=-log(d1p);
        aminoCont.dfy(1,:,:)=-log(d2p);
        
        %figure
        imagesc(sqrt(d1p.^2+d2p.^2));
        title([aminoName ':norm(Grad) phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
        
        
        print(gcf,'-dpsc2',strcat('images/Rama_Grad',aminoName,'.eps'))
        
        
        aminoCont.N_rota=col-7;
        
               
        save(['contData/' aminoName '.mat'],'aminoCont');
         
    end
    
        
end


