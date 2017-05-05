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

binAng = 360;

bin=pi/binAng*8.6;
kappa=5;

NumBins=round(2*pi/bin)+1;

bins=@(omega) round((pi+omega)/bin)+1;

%filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(abs(dC_CA-1.99)>1.18 || abs(dN_CA-1.78)>0.91 || abs(dP_plane-1.47)>0.55);
%filter=@(dC_CA,dN_CA,dP_plane,Bfactor) not(dC_CA<2*1.76 || dN_CA<2*1.7 || dP_plane<2*1.7);
%filter=@(dC_CA,dN_CA,dP_plane)1;
filter=@(dC_CA,dN_CA,dP_plane,Bfactor) 1;

%aminos={'PHE'; 'ASP'; 'THR'; 'ARG'; 'TRP'; 'VAL'; 'CYS'; 'SER'; 'ALA'; 'GLY'; 'MET'; 'TYR'; 'ASN'; 'PRO'; 'LYS'; 'HIS'; 'GLN'; 'ILE'; 'LEU'; 'GLU'};
aminos={'PHE'; 'ASP'}


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
    file = fopen(['_posts/2015-05-02-Amino-' aminoName  '.markdown'],'w');
    
    fprintf(file,'---\n');
    fprintf(file,'layout: post\n');
    fprintf(file,['title:  "Amino ' aminoName '"\n']);
    fprintf(file,['date:   2015-04-03 15:01:' num2str(i) ' +0100\n']);
    %fprintf(file,['date:   2015-04-03 15:01:' num2str(1) ' +0100\n']);
    fprintf(file,'categories: update\n');
    fprintf(file,'---\n');
    
    fprintf(file,'\n\n![Image](../../../../images/aadensity.png)\n');
    
    fprintf(file,['# Amino ' aminoName '\n']);
    
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
        
        fprintf(file,'\n\n # $$\\phi$$ & $$\\psi$$\n')
                
        figure
        [f,w1]=circ_vmpar(phi);
        [f,w2]=circ_vmpar(psi);
        [p nbins]=circ_vmpdf2(binAng, phi, psi, w1, w2);
        p=Rmatrix(p');
        imagesc(-log(p));
        title([aminoName ': psi vs phi']);
        xlabel('psi');
        ylabel('phi');
        %axis([-180 180 -180 180]);
        colorbar;
       
        print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_phipsi.jpg'));
        
        fprintf(file,['![Image](../../../../images/' aminoName '_Rama_phipsi.jpg)\n']);
                    
        figure
        
        [d1p,d2p] = imgradientxy(-log(p));
        imagesc(sqrt(d1p.^2+d2p.^2));
        title([aminoName ':norm(Grad) phi vs psi']);
        xlabel('psi');
        ylabel('phi');
        %axis([-180 180 -180 180]);
        colorbar;
             
        print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_phipsiGrad.jpg'))
        fprintf(file,['![Image](../../../../images/' aminoName '_Rama_phipsiGrad.jpg)\n']);
    
        
        
        if col> 8
                    chsi1= data(:,9);
                    
                    fprintf(file,'\n\n# $$\\chi_1$$\n')
                    
                    figure
                    [f,w3]=circ_vmpar(chsi1);
                    [p nbins]=circ_vmpdf2(binAng, phi, chsi1, w1,w3);
                    p=Rmatrix(p');
                    imagesc(-log(p));
                    title([aminoName ':phi vs chi1']);
                    xlabel('chi1');
                    ylabel('phi');
                    %axis([-180 180 -180 180]);
                    colorbar;
                    
                    print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_phichi1.jpg'))
                    fprintf(file,['![Image](../../../../images/' aminoName '_Rama_phichi1.jpg)\n']);
                    
        
                    figure
                    
                    [d1p,d2p] = imgradientxy(-log(p));
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) phi vs chi1']);
                    xlabel('chi1');
                    ylabel('phi');
                    %axis([-180 180 -180 180]);
                    colorbar;
                                     
                    print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_Grad_psichi1.jpg'))
                    fprintf(file,['![Image](../../../../images/' aminoName '_Rama_Grad_psichi1.jpg)\n']);

                    
                    figure
                    
                    [p nbins]=circ_vmpdf2(binAng, psi, chsi1, w2,w3);
                    %[p nbins]=circ_vmpdf2(binAng, psi, chsi1, w2,w3,p);
                    p=Rmatrix(p');
                    imagesc(-log(p));
                    title([aminoName ':psi vs chi1']);
                    xlabel('chi1');
                    ylabel('psi');
                    %axes('Position',[-180 180 -180 180]);
                    colorbar;
                    
                    print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_psichi1.jpg'))
                    fprintf(file,['![Image](../../../../images/' aminoName '_Rama_psichi1.jpg)\n']);
                    
        
                    figure

                    [d1p,d2p] = imgradientxy(-log(p));
                    aminoCont.dfx(3,:,:)=d1p;
                    aminoCont.dfy(3,:,:)=d2p;
                    imagesc(sqrt(d1p.^2+d2p.^2));
                    title([aminoName ':norm(Grad) psi vs chi1']);
                    xlabel('psi');
                    ylabel('chi1');
                    %axis([-180 180 -180 180]);
                    colorbar;
                    
                    
                    print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_Grad_phichi1.jpg'))
                    fprintf(file,['![Image](../../../../images/' aminoName '_Rama_Grad_phichi1.jpg)\n']);

                    
                    if col > 9                        
                        chsi2= data(:,10);
                        
                        fprintf(file,'\n\n# $$\\chi_2$$\n')
                        
                        figure
                        [f,w4]=circ_vmpar(chsi2);
                        [p nbins]=circ_vmpdf2(binAng, chsi1, chsi2, w3,w4);
                        p=Rmatrix(p');
                        imagesc(-log(p));
                        title([aminoName ':chi1 vs chi2']);
                        xlabel('chi2');
                        ylabel('chi1');
                        %axis([-180 180 -180 180]);
                        colorbar;    
                        
                        print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_chi1chi2.jpg'))
                        fprintf(file,['![Image](../../../../images/' aminoName '_Rama_chi1chi2.jpg)\n']);
                   

                        figure

                        [d1p,d2p] = imgradientxy(-log(p));
                        imagesc(sqrt(d1p.^2+d2p.^2));
                        title([aminoName ':norm(Grad) chi1 vs chi2']);
                        xlabel('chi2');
                        ylabel('chi1');
                        %axis([-180 180 -180 180]);
                        colorbar;
                        
                        print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_Gradchi1chi2.jpg'))
                        fprintf(file,['![Image](../../../../images/' aminoName '_Rama_Gradchi1chi2.jpg)\n']);
                    
                    
                        if col > 10
                            
                            chsi3= data(:,11);
                            
                            fprintf(file,'\n\n# $$\\chi_3$$\n')

                            figure
                            [f,w5]=circ_vmpar(chsi3);
                            [p nbins]=circ_vmpdf2(binAng, chsi2, chsi3, w4,w5);
                            p=Rmatrix(p');
                            imagesc(-log(p))
                            title([aminoName ':chi2 vs chi3']);
                            xlabel('chi3');
                            ylabel('chi2');
                            %axis([-180 180 -180 180]);
                            colorbar;
                            
                            print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_chsi2chi3.jpg'))
                            fprintf(file,['![Image](../../../../images/' aminoName '_Rama_chsi2chi3.jpg)\n']);
                    

                            figure

                            [d1p,d2p] = imgradientxy(-log(p));
                            imagesc(sqrt(d1p.^2+d2p.^2));
                            title([aminoName ':norm(Grad) chi2 vs chi3']);
                            xlabel('chi3');
                            ylabel('chi2');
                            %axis([-180 180 -180 180]);
                            colorbar;    
                            
                            print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_Gradchsi2chi3.jpg'))
                            fprintf(file,['![Image](../../../../images/' aminoName '_Rama_Gradchsi2chi3.jpg)\n']);
                    
                            
                            if col > 11
                                
                                figure
                                chsi4= data(:,12);

                                fprintf(file,'\n\n# $$\\chi_4$$\n')
                                [f,w6]=circ_vmpar(chsi4);
                                [p nbins]=circ_vmpdf2(binAng, chsi3, chsi4, w5,w6);
                                p=Rmatrix(p');
                                imagesc(-log(p));
                                title([aminoName ':chi3 vs chi4']);
                                xlabel('chi4');
                                ylabel('chi3');
                                %axis([-180 180 -180 180]);
                                colorbar;
                                
                                print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_chi3chi4.jpg'))
                                fprintf(file,['![Image](../../../../images/' aminoName '_Rama_chi3chi4.jpg)\n']);
                    


                                figure

                                [d1p,d2p] = imgradientxy(-log(p));
                                imagesc(sqrt(d1p.^2+d2p.^2));
                                title([aminoName ':norm(Grad) chi3 vs chi4']);
                                xlabel('chi4');
                                ylabel('chi3');
                                %axis([-180 180 -180 180]);
                                colorbar;
                                
                                print(gcf,'-djpeg',strcat('images/',aminoName,'_Rama_Gradchi3chi4.jpg'))
                                fprintf(file,['![Image](../../../../images/' aminoName '_Rama_Gradchi3chi4.jpg)\n']);
                    
                                
                            end
                        end
                    end
        end
        
        
%         %%
%         % Amino discrete characteristic
%         %
%         %% 
%         
%         data=[bins(Data(:,1)) bins(Data(:,2)) Data(:,3:end)];
%                 
%         amino.count=zeros(NumBins,NumBins,8);
%         
%         for i=1:NumBins
%             for j=1:NumBins
%                 infor.dC_CA=[];
%                 infor.dN_CA=[];
%                 infor.R1=[];
%                 infor.R2=[];
%                 infor.R3=[];
%                 infor.R4=[];
%                 amino=setfield(amino,['R' int2str(i) 'x' int2str(j)],infor);
%             end
%         end
%         
%         
%         for idx=1:l
%             i=data(idx,1);
%             j=data(idx,2);
%             k=1;
%             if ~(isnan(i) || isnan(j)) %&& filter(data(idx,3),data(idx,4),data(idx,6))
%                 amino.count(i,j,1)=amino.count(i,j,1)+1;
% 
%                 BIN = ['R' int2str(i) 'x' int2str(j)];
% 
%                 dataBIN = getfield(amino,BIN); 
% 
%                 dataBIN.dC_CA=[dataBIN.dC_CA data(idx,3)];
%                 dataBIN.dN_CA=[dataBIN.dN_CA data(idx,4)];
%                 if col>7 && ~isnan(data(idx,8))
%                     dataBIN.R1=[dataBIN.R1 data(idx,8)];
%                     k1=bins(data(idx,8));
%                     amino.count(i,k1,2)=amino.count(i,k1,2)+1;
%                     amino.count(k1,j,3)=amino.count(k1,j,3)+1;
%                     if col>8 && ~isnan(data(idx,9))
%                         dataBIN.R2=[dataBIN.R2 data(idx,9)];
%                         k2=bins(data(idx,9));
%                         amino.count(i,k2,4)=amino.count(i,k2,4)+1;
%                         amino.count(k2,j,5)=amino.count(k2,j,5)+1;
%                         amino.count(k1,k2,6)=amino.count(k1,k2,6)+1;
%                         if col>9 && ~isnan(data(idx,10))
%                             dataBIN.R3=[dataBIN.R3 data(idx,10)];
%                             k3=bins(data(idx,10));
%                             amino.count(k1,k3,7)=amino.count(k1,k3,7)+1;
%                             amino.count(k2,k3,8)=amino.count(k2,k3,8)+1;
%                             if col> 10 && ~isnan(data(idx,11))
%                                 dataBIN.R4=[dataBIN.R4 data(idx,11)];
%                             end
%                         end
%                     end
%                 end
%                 amino=setfield(amino,BIN,dataBIN);
%             end
%         end
%         
%         amino=setfield(amino,'Size',l);
%         amino=setfield(amino,'NumBins',NumBins);
%         
%         amino.MeandC_CA=zeros(NumBins,NumBins);
%         amino.VardC_CA=zeros(NumBins,NumBins);
%         amino.mindC_CA=zeros(NumBins,NumBins);
%         amino.MeanddN_CA=zeros(NumBins,NumBins);
%         amino.VardN_CA=zeros(NumBins,NumBins);
%         amino.mindN_CA=zeros(NumBins,NumBins);
%         if col>7
%             amino.MeanR1=zeros(NumBins,NumBins);
%             amino.VarR1=zeros(NumBins,NumBins);
%             if col>8
%                 amino.MeanR2=zeros(NumBins,NumBins);
%                 amino.VarR2=zeros(NumBins,NumBins);
%                 if col>9
%                     amino.MeanR3=zeros(NumBins,NumBins);
%                     amino.VarR3=zeros(NumBins,NumBins);
%                     if col> 10
%                         amino.MeanR4=zeros(NumBins,NumBins);
%                         amino.VarR4=zeros(NumBins,NumBins);
%                     end
%                 end
%             end
%         end
%         for i=1:NumBins
%             for j=1:NumBins
%                 BIN = ['R' int2str(i) 'x' int2str(j)];            
%                 dataBIN = getfield(amino,BIN);
%                 
%                 amino.MeandC_CA(i,j)=mean(dataBIN.dC_CA);
%                 amino.VardC_CA(i,j)=var(dataBIN.dC_CA);
%                 
%                 m=min(dataBIN.dC_CA);
%                 if ~isempty(m)
%                     amino.minC_CA(i,j)=m;
%                 else
%                     amino.minC_CA(i,j)=nan;
%                 end
%                 
%                 amino.MeandN_CA(i,j)=mean(dataBIN.dN_CA);
%                 amino.VardN_CA(i,j)=var(dataBIN.dN_CA);
%                 
%                 m=min(dataBIN.dN_CA);
%                 if ~isempty(m)
%                     amino.minN_CA(i,j)=m;
%                 else
%                     amino.minN_CA(i,j)=nan;
%                 end
%                 
%                 if col>7
%                     amino.MeanR1(i,j)=mean(dataBIN.R1);
%                     amino.VarR1(i,j)=var(dataBIN.R1);
%                     if col>8
%                         amino.MeanR2(i,j)=mean(dataBIN.R2);
%                         amino.VarR2(i,j)=var(dataBIN.R2);
%                         if col>9
%                             amino.MeanR3(i,j)=mean(dataBIN.R3);
%                             amino.VarR3(i,j)=var(dataBIN.R3);
%                             if col> 10
%                                 amino.MeanR4(i,j)=mean(dataBIN.R4);
%                                 amino.VarR4(i,j)=var(dataBIN.R4);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
    end
    
        fclose('all');
end


