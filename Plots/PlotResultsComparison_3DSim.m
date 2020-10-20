function [] = PlotResultsComparison_3DSim(ResultsFile1,ResultsFile2)

Csv = hsv(3);

h = figure();
h.Name = 'Sim3D_ResultsComparison';
hTabGroup = uitabgroup;
tab1 = uitab(hTabGroup, 'Title', 'Ankle 1');
tab2 = uitab(hTabGroup, 'Title', 'Ankle 2');
set(h,'Color','w');

            
for j=1:3
if(j==1)
    [path,file,~] = fileparts(ResultsFile1);
    addpath(path);
    load(file);
elseif(j==2)
   [path,file,~] = fileparts(ResultsFile2);
    addpath(path);
    load(file);
end

Cs = Csv(j,:);

    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
    iTib = find(strcmp(R.colheaders.muscles,'tib_ant_r'));
    if isempty(iGas)
        iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
    end
    if isempty(iGas2)
        iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
    end
    if isempty(iTib)
        iTib = find(strcmp(R.colheaders.muscles,'tibant_r'));
    end
    mVect = {'Soleus','Gas-lat','Gas-med','Tib-ant'};
    
    iM = [iSol iGas iGas2 iTib];
    
    for i=1:4
        if(j==1)
            a1 = R.a;
            MetabB_Etot1 = R.MetabB.Etot;
            lMtilde1 = R.lMtilde;
            FT1 = R.FT;
        end
        
        if(j==3)    % plot differences
            axes('parent', tab1);
        subplot(4,4,i+4)
        plot((R.a(:,iM(i))-a1(:,iM(i))),'-','Color',Cs); hold on; grid on
        xlabel('% stride'); ylabel('\Delta activity');
        
        subplot(4,4,i+12)
        plot((R.MetabB.Etot(:,iM(i))-MetabB_Etot1(:,iM(i))),'-','Color',Cs); hold on; grid on;
        xlabel('% stride'); ylabel('\Delta Muscle metab power');
        
            axes('parent', tab2);
        subplot(4,4,i+4)
        plot((R.lMtilde(:,iM(i))-lMtilde1(:,iM(i))),'-','Color',Cs); hold on; grid on;
        xlabel('% stride'); ylabel('\Delta Norm fiber length');
        
        subplot(4,4,i+12)
        plot((R.FT(:,iM(i))-FT1(:,iM(i))),'-','Color',Cs); hold on; grid on;
        xlabel('% stride'); ylabel('\Delta Norm muscle force');
        
        else    % plot values
            axes('parent', tab1);
        subplot(4,4,i)
        plot(R.a(:,iM(i)),'-','Color',Cs); hold on;  title(mVect{i});
        xlabel('% stride'); ylabel('activity');
        
        subplot(4,4,i+8)
        plot(R.MetabB.Etot(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Muscle metab power');
        
            axes('parent', tab2);
        subplot(4,4,i)
        plot(R.lMtilde(:,iM(i)),'-','Color',Cs); hold on; title(mVect{i});
        xlabel('% stride'); ylabel('Norm fiber length');
        
        subplot(4,4,i+8)
        plot(R.FT(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Norm muscle force');
        end
    end
    
   
    
end
    
    
    