%%
% 
% Analysis: A simple script used to scan the extracted data for each
% aminoacid ploting distribution graphs for dihedral angles and bond
% lengths
%
%  Dependences: importdata, removenan, ksdensity, circ_vmpdf2
%  libs:        pdbTools, CircStat

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

filter=@(dC_CA,dN_CA,dP_plane) not(abs(dC_CA-1.99)>1.18 || abs(dN_CA-1.78)>0.91 || abs(dP_plane-1.47)>0.55);
%filter=@(dC_CA,dN_CA,dP_plane) not(dC_CA<2*1.76 || dN_CA<2*1.7 || dP_plane<2*1.7);
%filter=@(dC_CA,dN_CA,dP_plane)1;

aminos={'PHE'; 'ASP'; 'THR'; 'ARG'; 'TRP'; 'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
%aminos={'ARG'}

% Attributes extrated from the protein db:
% dihedral angles............. phi, psi
% bond distance............... dC_CA, dN_CA, dN_C
% plolypetid plane............ dP_plan
% angle....................... ang_N_CA_C
% side chain dihedral angles.. R1, R2, R3, R4
varName={'phi', 'psi', 'dC_CA','dN_CA','dN_C','dP_plane', 'ang_N_CA_C','R1','R2','R3','R4'};

for i=1:length(aminos)   
    aminoName=aminos{i}
    % import the data in CSV format used to the statistical analysies
    %
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
                
        figure
        [f,xi,w1]=ksdensity(phi);
        [f,xi,w2]=ksdensity(psi);
        kappa=min([w1,w2]);
        %lambda = circ_density( phi,psi, kappa)
        [p nbins]=circ_vmpdf2(180, phi, psi, kappa);
        imagesc(-log(p));
        title([aminoName ': phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
       
        print(gcf,'-dpsc2',strcat('images/Rama_',aminoName,'.eps'))
                    
        
        figure
        [d1 nbins] = circ_d1vmpdf2(180, phi,psi, kappa);
        [d2 nbins] = circ_d2vmpdf2(180, phi,psi, kappa);
        d1p=-d1./p;
        d2p=-d2./p;
        imagesc(sqrt(d1p.^2+d2p.^2));
        title([aminoName ':norm(Grad) phi vs psi']);
        xlabel('phi');
        ylabel('psi');
        %%axis([-180 180 -180 180]);
        colorbar;
        
        
        print(gcf,'-dpsc2',strcat('images/Rama_Grad',aminoName,'.eps'))
        
        %aminoCont.dphi_psi=d1p:            
        %aminoCont.phi_dpsi=d2p:
        
        
        if col> 7
                    chsi1= data(:,8);
                    
                    fprintf('\t\t\tCHI_1\n')
                    
                    figure
                    [f,xi,w3]=ksdensity(chsi1);
                    [p nbins]=circ_vmpdf2(180, phi, chsi1, (w1+w3)/2);
                    imagesc(-log(p));
                    title([aminoName ':phi vs chi1']);
                    xlabel('phi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    print(gcf,'-dpsc2',strcat('images/Rama_phi_chi1',aminoName,'.eps'))
                    
        
                    figure
                    
                    [d1 nbins] = circ_d1vmpdf2(180, phi, chsi1, (w1+w3)/2);
                    [d2 nbins] = circ_d2vmpdf2(180, phi, chsi1, (w1+w3)/2);
                    d1p=-d1./p;
                    d2p=-d2./p;
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) phi vs chi1']);
                    xlabel('phi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    
                    print(gcf,'-dpsc2',strcat('images/Rama_Grad_psichi1',aminoName,'.eps'))

                    figure
                    [p nbins]=circ_vmpdf2(180, psi, chsi1, (w2+w3)/2);
                    imagesc(-log(p));
                    title([aminoName ':psi vs chi1']);
                    xlabel('psi');
                    ylabel('chi1');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    print(gcf,'-dpsc2',strcat('images/Rama_psi_chi1',aminoName,'.eps'))
                    
        
                    figure
                    [d1 nbins] = circ_d1vmpdf2(180, psi, chsi1, (w2+w3)/2);
                    [d2 nbins] = circ_d2vmpdf2(180, psi, chsi1, (w2+w3)/2);
                    d1p=-d1./p;
                    d2p=-d2./p;
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) psi vs chi1']);
                    xlabel('chi1');
                    ylabel('psi');
                    %%axis([-180 180 -180 180]);
                    colorbar;
                    
                    
                    print(gcf,'-dpsc2',strcat('images/Rama_Grad_phichi1',aminoName,'.eps'))

                    
                    if col > 8                        
                        chsi2= data(:,9);
                        
                        fprintf('\t\t\tCHI_2\n')
                        
                        figure
                        [f,xi,w4]=ksdensity(chsi2);
                        [p nbins]=circ_vmpdf2(180, chsi1, chsi2, (w4+w3)/2);
                        imagesc(-log(p));
                        title([aminoName ':chi1 vs chi2']);
                        xlabel('chi1');
                        ylabel('chi2');
                        %%axis([-180 180 -180 180]);
                        colorbar;    
                        
                        print(gcf,'-dpsc2',strcat('images/Rama_chi2',aminoName,'.eps'))
                   

                        figure
                        [d1 nbins] = circ_d1vmpdf2(180, chsi1, chsi2, (w4+w3)/2);
                        [d2 nbins] = circ_d2vmpdf2(180, chsi1, chsi2, (w4+w3)/2);
                        d1p=-d1./p;
                        d2p=-d2./p;
                        imagesc(sqrt(d1p.^2+d2p.^2));
                        title([aminoName ':norm(Grad) chi1 vs chi2']);
                        xlabel('chi1');
                        ylabel('chi2');
                        %%axis([-180 180 -180 180]);
                        colorbar;
                        
                        print(gcf,'-dpsc2',strcat('images/Rama_Gradchi2',aminoName,'.eps'))
                    
                    
                        if col > 9
                            
                            chsi3= data(:,10);
                            
                            fprintf('\t\t\tCHI_3\n')

                            figure
                            [f,xi,w5]=ksdensity(chsi3);
                            [p nbins]=circ_vmpdf2(180, chsi2, chsi3, (w4+w5)/2);
                            imagesc(-log(p))
                            title([aminoName ':chi2 vs chi3']);
                            xlabel('chi2');
                            ylabel('chi3');
                            %%axis([-180 180 -180 180]);
                            colorbar;
                            
                            print(gcf,'-dpsc2',strcat('images/Rama_chi3',aminoName,'.eps'))
                    

                            figure
                            [d1 nbins] = circ_d1vmpdf2(180, chsi2, chsi3, (w4+w5)/2);
                            [d2 nbins] = circ_d2vmpdf2(180, chsi2, chsi3, (w4+w5)/2);
                            d1p=-d1./p;
                            d2p=-d2./p;
                            imagesc(sqrt(d1p.^2+d2p.^2));
                            title([aminoName ':norm(Grad) chi2 vs chi3']);
                            xlabel('chi2');
                            ylabel('chi3');
                            %%axis([-180 180 -180 180]);
                            colorbar;    
                            
                            print(gcf,'-dpsc2',strcat('images/Rama_Gradchi3',aminoName,'.eps'))
                    
                            
                            if col > 10
                                
                                chsi4= data(:,11);
                                
                                fprintf('\t\t\tCHI_4\n')
                                
                                figure
                                [f,xi,w6]=ksdensity(chsi4);
                                [p nbins]=circ_vmpdf2(180, chsi3, chsi4, kappa);
                                imagesc(-log(p));
                                title([aminoName ':chi3 vs chi4']);
                                xlabel('chi3');
                                ylabel('chi4');
                                %%axis([-180 180 -180 180]);
                                colorbar;
                                
                                print(gcf,'-dpsc2',strcat('images/Rama_chi4',aminoName,'.eps'))
                    


                                figure
                                [d1 nbins] = circ_d1vmpdf2(180, chsi3, chsi4, (w6+w5)/2);
                                [d2 nbins] = circ_d2vmpdf2(180, chsi3, chsi4, (w6+w5)/2);
                                d1p=-d1./p;
                                d2p=-d2./p;
                                imagesc(sqrt(d1p.^2+d2p.^2));
                                title([aminoName ':norm(Grad) chi3 vs chi4']);
                                xlabel('chi3');
                                ylabel('chi4');
                                %%axis([-180 180 -180 180]);
                                colorbar;
                                
                                print(gcf,'-dpsc2',strcat('images/Rama_Gradchi4',aminoName,'.eps'))
                    
                                
                            end
                        end
                    end
        end
        
        
        save(['binData/' aminoName '.mat'],'aminoCont');
        
        
        %%
        % Amino discrete characteristic
        %
        
        data=[bins(Data(:,1)) bins(Data(:,2)) Data(:,3:end)];
                
        amino.count=zeros(NumBins,NumBins,8);
        
        for i=1:NumBins
            for j=1:NumBins
                infor.dC_CA=[];
                infor.dN_CA=[];
                infor.R1=[];
                infor.R2=[];
                infor.R3=[];
                infor.R4=[];
                amino=setfield(amino,['R' int2str(i) 'x' int2str(j)],infor);
            end
        end
        
        
        for idx=1:l
            i=data(idx,1);
            j=data(idx,2);
            k=1;
            if ~(isnan(i) || isnan(j)) %&& filter(data(idx,3),data(idx,4),data(idx,6))
                amino.count(i,j,1)=amino.count(i,j,1)+1;

                BIN = ['R' int2str(i) 'x' int2str(j)];

                dataBIN = getfield(amino,BIN); 

                dataBIN.dC_CA=[dataBIN.dC_CA data(idx,3)];
                dataBIN.dN_CA=[dataBIN.dN_CA data(idx,4)];
                if col>7 && ~isnan(data(idx,8))
                    dataBIN.R1=[dataBIN.R1 data(idx,8)];
                    k1=bins(data(idx,8));
                    amino.count(i,k1,2)=amino.count(i,k1,2)+1;
                    amino.count(k1,j,3)=amino.count(k1,j,3)+1;
                    if col>8 && ~isnan(data(idx,9))
                        dataBIN.R2=[dataBIN.R2 data(idx,9)];
                        k2=bins(data(idx,9));
                        amino.count(i,k2,4)=amino.count(i,k2,4)+1;
                        amino.count(k2,j,5)=amino.count(k2,j,5)+1;
                        amino.count(k1,k2,6)=amino.count(k1,k2,6)+1;
                        if col>9 && ~isnan(data(idx,10))
                            dataBIN.R3=[dataBIN.R3 data(idx,10)];
                            k3=bins(data(idx,10));
                            amino.count(k1,k3,7)=amino.count(k1,k3,7)+1;
                            amino.count(k2,k3,8)=amino.count(k2,k3,8)+1;
                            if col> 10 && ~isnan(data(idx,11))
                                dataBIN.R4=[dataBIN.R4 data(idx,11)];
                            end
                        end
                    end
                end
                amino=setfield(amino,BIN,dataBIN);
            end
        end
        
        amino=setfield(amino,'Size',l);
        amino=setfield(amino,'NumBins',NumBins);
        
        amino.MeandC_CA=zeros(NumBins,NumBins);
        amino.VardC_CA=zeros(NumBins,NumBins);
        amino.mindC_CA=zeros(NumBins,NumBins);
        amino.MeanddN_CA=zeros(NumBins,NumBins);
        amino.VardN_CA=zeros(NumBins,NumBins);
        amino.mindN_CA=zeros(NumBins,NumBins);
        if col>7
            amino.MeanR1=zeros(NumBins,NumBins);
            amino.VarR1=zeros(NumBins,NumBins);
            if col>8
                amino.MeanR2=zeros(NumBins,NumBins);
                amino.VarR2=zeros(NumBins,NumBins);
                if col>9
                    amino.MeanR3=zeros(NumBins,NumBins);
                    amino.VarR3=zeros(NumBins,NumBins);
                    if col> 10
                        amino.MeanR4=zeros(NumBins,NumBins);
                        amino.VarR4=zeros(NumBins,NumBins);
                    end
                end
            end
        end
        for i=1:NumBins
            for j=1:NumBins
                BIN = ['R' int2str(i) 'x' int2str(j)];            
                dataBIN = getfield(amino,BIN);
                
                amino.MeandC_CA(i,j)=mean(dataBIN.dC_CA);
                amino.VardC_CA(i,j)=var(dataBIN.dC_CA);
                
                m=min(dataBIN.dC_CA);
                if ~isempty(m)
                    amino.minC_CA(i,j)=m;
                else
                    amino.minC_CA(i,j)=nan;
                end
                
                amino.MeandN_CA(i,j)=mean(dataBIN.dN_CA);
                amino.VardN_CA(i,j)=var(dataBIN.dN_CA);
                
                m=min(dataBIN.dN_CA);
                if ~isempty(m)
                    amino.minN_CA(i,j)=m;
                else
                    amino.minN_CA(i,j)=nan;
                end
                
                if col>7
                    amino.MeanR1(i,j)=mean(dataBIN.R1);
                    amino.VarR1(i,j)=var(dataBIN.R1);
                    if col>8
                        amino.MeanR2(i,j)=mean(dataBIN.R2);
                        amino.VarR2(i,j)=var(dataBIN.R2);
                        if col>9
                            amino.MeanR3(i,j)=mean(dataBIN.R3);
                            amino.VarR3(i,j)=var(dataBIN.R3);
                            if col> 10
                                amino.MeanR4(i,j)=mean(dataBIN.R4);
                                amino.VarR4(i,j)=var(dataBIN.R4);
                            end
                        end
                    end
                end
            end
        end
        
        amino.max_bin=max(max(amino.count(:,:,1)));
        
        amino.String=1/2*KC_C*(amino.MeandC_CA-MeandC_CA)^2 + 1/2*KN_C*(amino.MeandN_CA-MeandN_CA)^2;
        
        
        save(['binData/' aminoName '.mat'],'amino');
        
        %imagesc(amino.count(:,1,:));
        
        %title(['Ramachandran:' aminoName]);
        %xlabel('phi');
        %ylabel('psi');
        
        %print(gcf,'-dpsc2',strcat('images/Rama_',aminoName,'.eps'))
        
        %pause
        
        %imagesc(-log(1+amino.count(:,:,2)));
        
        %title(['Free energy:' aminoName]);
        %xlabel('phi');
        %ylabel('psi');
        
        %print(gcf,'-dpsc2',strcat('images/FreeEner_',aminoName,'.eps'))
        
        %pause
        figure
        imagesc(amino.count(:,:,1)/amino.max_bin);
        
        figure
        [p nbins]=circ_vmpdf2(45, phi, psi, 1);
        imagesc(p);
        
        pause
        imagesc(amino.count(:,:,2)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,3)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,4)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,5)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,6)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,7)/amino.max_bin);
        pause
        imagesc(amino.count(:,:,8)/amino.max_bin);
        pause
        
    end
    
        
end


