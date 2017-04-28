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
aminos={'ARG'}
varName={'phi', 'psi', 'dC_CA','dN_CA','dN_C','dP_plane', 'ang_N_CA_C','R1','R2','R3','R4'};

for i=1:length(aminos)   
    aminoName=aminos{i}
    
    Data=importdata(['aminoData/' aminoName '.cvs']);
    
    if ~isempty(Data)
        
        [l,col]=size(Data);
        
        data=removenan(Data,filter); % clean and select data
        
        %%
        % Amino discrete characteristic
        %
        phi= data(:,1);
        psi= data(:,2);
        
        % bounds
        C_CA= data(:,3);
        N_CA= data(:,4);
        dplane= data(:,6);
                
        
        [p nbins]=circ_harmonic(180,phi, psi);
        %figure
        imagesc(exp(-p));
        
        %pause
        
        aminoCont.f(1,:,:)=exp(-p);
        
        title([aminoName ': phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
       
        %print(gcf,'-dpsc2',strcat('images/Rama_',aminoName,'.eps'))
                    
        
        [d1 nbins] = circ_d1harmonic(180,phi, psi);
        [d2 nbins] = circ_d2harmonic(180,phi, psi);
        aminoCont.dfx(1,:,:)=exp(-d1);
        aminoCont.dfy(1,:,:)=exp(-d2);
        
        %figure
        imagesc(sqrt(d1.^2+d2.^2));
        title([aminoName ':norm(Grad) phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
        
        
        %print(gcf,'-dpsc2',strcat('images/Rama_Grad',aminoName,'.eps'))
        
        
        aminoCont.N_rota=col-7;
        
        if col> 7
                    chsi1= data(:,8);
                    
                    fprintf('\t\t\tCHI_1\n')
                    
                    [f,w3]=circ_vmpar(chsi1);
                    [p nbins]=circ_harmonic(180,phi, chsi1); 
                    %figure
                    imagesc(exp(-p));
                    aminoCont.f(2,:,:)=exp(-p);
                    title([aminoName ':phi vs chi1']);
                    xlabel('phi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    %print(gcf,'-dpsc2',strcat('images/Rama_phi_chi1',aminoName,'.eps'))
                    
        
                    [d1 nbins] = circ_d1harmonic(180,phi, chsi1);
                    [d2 nbins] = circ_d2harmonic(180,phi, chsi1);
                    aminoCont.dfx(2,:,:)=exp(-d1p);
                    aminoCont.dfy(2,:,:)=exp(-d2p);
        
                    %figure
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) phi vs chi1']);
                    xlabel('phi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    
                    %print(gcf,'-dpsc2',strcat('images/Rama_Grad_psichi1',aminoName,'.eps'))

                    [p nbins]=circ_harmonic(180,psi, chsi1); 
                    %figure
                    imagesc(exp(-p));
                    aminoCont.f(3,:,:)=exp(-p);
                    title([aminoName ':psi vs chi1']);
                    xlabel('psi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    %print(gcf,'-dpsc2',strcat('images/Rama_psi_chi1',aminoName,'.eps'))
                    

                    [d1 nbins] = circ_d1harmonic(180,psi, chsi1); 
                    [d2 nbins] = circ_d2harmonic(180,psi, chsi1);
                    aminoCont.dfx(3,:,:)=exp(-d1p);
                    aminoCont.dfy(3,:,:)=exp(-d2p);
                    
                    %figure
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) psi vs chi1']);
                    xlabel('chi1');
                    ylabel('psi');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    
                    %print(gcf,'-dpsc2',strcat('images/Rama_Grad_phichi1',aminoName,'.eps'))

                    
                    if col > 8                        
                        chsi2= data(:,9);
                        
                        fprintf('\t\t\tCHI_2\n')
                        
                        [f,w4]=circ_vmpar(chsi2);
                        [p nbins]=circ_harmonic(180,chsi1, chsi2);
                        %figure
                        imagesc(exp(-p));
                        aminoCont.f(4,:,:)=exp(-p);
                        title([aminoName ':chi1 vs chi2']);
                        xlabel('chi1');
                        ylabel('chi2');
                        %%axis([-180 180 -180 180]);
                        colorbar;    
                        
                        %print(gcf,'-dpsc2',strcat('images/Rama_chi2',aminoName,'.eps'))
                   

                        %figure
                        [d1 nbins] = circ_d1harmonic(180,chsi1, chsi2);
                        [d2 nbins] = circ_d2harmonic(180,chsi1, chsi2);
                        aminoCont.dfx(4,:,:)=exp(-d1);
                        aminoCont.dfy(4,:,:)=exp(-d2);
                        
                        imagesc(sqrt(d1.^2+d2.^2));
                        title([aminoName ':norm(Grad) chi1 vs chi2']);
                        xlabel('chi1');
                        ylabel('chi2');
                        %%axis([-180 180 -180 180]);
                        colorbar;
                        
                        %print(gcf,'-dpsc2',strcat('images/Rama_Gradchi2',aminoName,'.eps'))
                    
                    
                        if col > 9
                            
                            chsi3= data(:,10);
                            
                            fprintf('\t\t\tCHI_3\n')

                            %figure
                            [f,w5]=circ_vmpar(chsi3);
                            [p nbins]=circ_harmonic(180,chsi2, chsi3);
                            imagesc(exp(-p))
                            aminoCont.f(5,:,:)=exp(-p);
                            title([aminoName ':chi2 vs chi3']);
                            xlabel('chi2');
                            ylabel('chi3');
                            %%axis([-180 180 -180 180]);
                            colorbar;
                            
                            %print(gcf,'-dpsc2',strcat('images/Rama_chi3',aminoName,'.eps'))
                    

                            %figure
                            [d1 nbins] = circ_d1harmonic(180,chsi2, chsi3);
                            [d2 nbins] = circ_d2harmonic(180,chsi2, chsi3);
                            aminoCont.dfx(5,:,:)=exp(-d1);
                            aminoCont.dfy(5,:,:)=exp(-d2);
                            
                            imagesc(sqrt(d1.^2+d2.^2));
                            title([aminoName ':norm(Grad) chi2 vs chi3']);
                            xlabel('chi2');
                            ylabel('chi3');
                            %%axis([-180 180 -180 180]);
                            colorbar;    
                            
                            %print(gcf,'-dpsc2',strcat('images/Rama_Gradchi3',aminoName,'.eps'))
                    
                            
                            if col > 10
                                
                                chsi4= data(:,11);
                                
                                fprintf('\t\t\tCHI_4\n')
                                
                                %figure
                                [f,w6]=circ_vmpar(chsi4);
                                [p nbins]=circ_harmonic(180,chsi3, chsi4); 
                                imagesc(exp(-p));
                                aminoCont.f(6,:,:)=exp(-p);
                                title([aminoName ':chi3 vs chi4']);
                                xlabel('chi3');
                                ylabel('chi4');
                                %%axis([-180 180 -180 180]);
                                colorbar;
                                
                                %print(gcf,'-dpsc2',strcat('images/Rama_chi4',aminoName,'.eps'))
                    


                                %figure
                                [d1 nbins] = circ_d1harmonic(180,chsi3, chsi4);
                                [d2 nbins] = circ_d2harmonic(180,chsi3, chsi4);
                                aminoCont.dfx(6,:,:)=exp(-d1);
                                aminoCont.dfy(6,:,:)=exp(-d2);
                                
                                imagesc(sqrt(d1.^2+d2.^2));
                                title([aminoName ':norm(Grad) chi3 vs chi4']);
                                xlabel('chi3');
                                ylabel('chi4');
                                %%axis([-180 180 -180 180]);
                                colorbar;
                                
                                %print(gcf,'-dpsc2',strcat('images/Rama_Gradchi4',aminoName,'.eps'))
                    
                                
                            end
                        end
                    end
        end
        
        
        %save(['contData/' aminoName '.mat'],'aminoCont');
         
    end
    
        
end


