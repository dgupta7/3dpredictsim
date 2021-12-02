function [] = PlotEnergySmoothing(ResultsFile,LegNames)

nr = length(ResultsFile);
if numel(LegNames) == nr
    LN = 1;
else
    LN = 0;
end
set(0,'defaultTextInterpreter','tex');

CsV = hsv(nr);
scs = get(0,'ScreenSize');
fsc = scs;
fsc(2) = 20;
fsc(4) = fsc(4)-120;

h = figure('Position',fsc);
h.Name = 'Energy model smoothing';
hTabGroup = uitabgroup;
tab1 = uitab(hTabGroup, 'Title', '1');
tab2 = uitab(hTabGroup, 'Title', '2');
tab3 = uitab(hTabGroup, 'Title', '3');
tab4 = uitab(hTabGroup, 'Title', '4');
tab5 = uitab(hTabGroup, 'Title', '5');
tab6 = uitab(hTabGroup, 'Title', '6');
tab7 = uitab(hTabGroup, 'Title', '7');
tab8 = uitab(hTabGroup, 'Title', '8');
tab9 = uitab(hTabGroup, 'Title', '9');
tab10 = uitab(hTabGroup, 'Title', '10');
tab11 = uitab(hTabGroup, 'Title', '11');
tab12 = uitab(hTabGroup, 'Title', '12');

iM = [47:92];


for inr=1:nr
    load(ResultsFile{inr},'R');
    if LN
        LegName = LegNames{inr};
    else
        LegName = R.S.savename;
    end
    Cs = CsV(inr,:);
   
    
    for i=1:size(iM,2)
        ii = rem(i-1,4)+1;
        if(i<5)
            axes('parent', tab1);
        elseif(i<9)
            axes('parent', tab2);    
        elseif(i<13)
            axes('parent', tab3);    
        elseif(i<17)
            axes('parent', tab4);
        elseif(i<21)
            axes('parent', tab5);
        elseif(i<25)
            axes('parent', tab6);
        elseif(i<29)
            axes('parent', tab7);
        elseif(i<33)
            axes('parent', tab8);    
        elseif(i<37)
            axes('parent', tab9);    
        elseif(i<41)
            axes('parent', tab10);    
        elseif(i<45)
            axes('parent', tab11);    
        else
            axes('parent', tab12);    
        end
            
%         subplot(3,4,ii)
%         plot(R.MetabB.Etot(:,iM(i)),'-','Color',Cs,'DisplayName',LegName); hold on;  
% %         title(mVect{i});
%         title(R.colheaders.muscles{iM(i)},'Interpreter','none')
% 
%         ylabel('Edot');
%         if ii==4
%             legend
%         end
%         
%         subplot(3,4,ii+4)
%         Hdot = R.MetabB.Adot(:,iM(i)) + R.MetabB.Mdot(:,iM(i)) + R.MetabB.Sdot(:,iM(i));
%         plot(Hdot,'-','Color',Cs); hold on;
%         ylabel('Hdot');
% 
%         subplot(3,4,ii+8)
%         plot(R.MetabB.Wdot(:,iM(i)),'-','Color',Cs); hold on;
%         Wdot = -R.Muscle.Fce(:,iM(i)).*R.Muscle.vM(:,iM(i));
%         plot(Wdot,':','Color',Cs)
%         ylabel('Wdot');
%         grid on
%         
%         xlabel('% stride')

        subplot(5,4,ii)
        plot(R.MetabB.Etot(:,iM(i)),'Color',Cs,'DisplayName',['Smoothed (b=' num2str(R.S.tanh_b) ')']); hold on;
        plot(R.MetabB_non_smooth.Etot(:,iM(i)),'--','Color',Cs,'DisplayName',['Non-smoothed (b=' num2str(R.S.tanh_b) ')'])
        title(R.colheaders.muscles{iM(i)},'Interpreter','none')
        ylabel('Edot'); grid on
        if ii==4
            lh = legend('Location','northwest');
            if inr==nr
                lhPos = lh.Position;
                lhPos(1) = lhPos(1)+0.15;
%                 lhPos(2) = lhPos(2)+0.1;
                set(lh,'position',lhPos);
            end
        end

        subplot(5,4,ii+4)
        plot(R.MetabB.Adot(:,iM(i)),'Color',Cs); hold on;
        plot(R.MetabB_non_smooth.Adot(:,iM(i)),'--','Color',Cs);
        ylabel('Adot'); grid on
    
        subplot(5,4,ii+8)
        plot(R.MetabB.Mdot(:,iM(i)),'Color',Cs); hold on;
        plot(R.MetabB_non_smooth.Mdot(:,iM(i)),'--','Color',Cs);
        ylabel('Mdot'); grid on
    
        subplot(5,4,ii+12)
        plot(R.MetabB.Sdot(:,iM(i)),'Color',Cs); hold on;
        plot(R.MetabB_non_smooth.Sdot(:,iM(i)),'--','Color',Cs);
        ylabel('Sdot'); grid on
    
        subplot(5,4,ii+16)
        plot(R.MetabB.Wdot(:,iM(i)),'Color',Cs); hold on;
        plot(R.MetabB_non_smooth.Wdot(:,iM(i)),'--','Color',Cs);
        ylabel('Wdot'); grid on

    end

    dist_trav = R.Qs(end,strcmp(R.colheaders.joints,'pelvis_tx')) - ...
                R.Qs(1,strcmp(R.colheaders.joints,'pelvis_tx'));


    getCOT  =@(Ed) trapz(R.t,sum(Ed,2))/R.body_mass/dist_trav;

    COT_E = getCOT(R.MetabB.Etot);
    COT_E_nsm = getCOT(R.MetabB_non_smooth.Etot);

    COT_A = getCOT(R.MetabB.Adot);
    COT_A_nsm = getCOT(R.MetabB_non_smooth.Adot);

    COT_M = getCOT(R.MetabB.Mdot);
    COT_M_nsm = getCOT(R.MetabB_non_smooth.Mdot);

    COT_S = getCOT(R.MetabB.Sdot);
    COT_S_nsm = getCOT(R.MetabB_non_smooth.Sdot);

    COT_W = getCOT(R.MetabB.Wdot);
    COT_W_nsm = getCOT(R.MetabB_non_smooth.Wdot);
    
    subplot(5,4,3)
    bar(inr*2-1,COT_E,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,COT_E_nsm,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('E (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';


    subplot(5,4,3+4)
    bar(inr*2-1,COT_A,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,COT_A_nsm,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('A (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
%     ylim([0,0.8])

    subplot(5,4,3+8)
    bar(inr*2-1,COT_M,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,COT_M_nsm,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('M (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
%     ylim([0,1.2])

    subplot(5,4,3+12)
    bar(inr*2-1,COT_S,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,COT_S_nsm,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('S (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
%     ylim([0,0.5])

    subplot(5,4,3+16)
    bar(inr*2-1,COT_W,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,COT_W_nsm,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('W (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
%     ylim([0,1.6])

    subplot(5,4,[4,8])
    bar(inr*2-1,R.COT,'FaceColor',Cs,'LineStyle','-'); hold on;
    bar(inr*2,R.COT_non_smooth,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5);
    ylabel('COT (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
%     ylim([0,5])

    subplot(5,4,[12,16])
    bar(inr*2-1,R.COT,'FaceColor',Cs,'LineStyle','-','DisplayName',['Smoothed (b=' num2str(R.S.tanh_b) ')']); hold on;
    bar(inr*2,R.COT_non_smooth,'FaceColor',Cs,'LineStyle','--','FaceAlpha',0.5,'DisplayName',['Non-smoothed (b=' num2str(R.S.tanh_b) ')']);
    ylabel('COT (Wkg^{-1}m^{-1})'); grid on
    ax1 = gca;
    ax1.XAxis.Visible = 'off';
    ylim([3.7,4.8])

    lh1 = legend('Location','northwest');
    if inr==nr
        lhPos = lh1.Position;
        lhPos(2) = lhPos(2)-0.3;
        set(lh1,'position',lhPos);
    end

end
    
end
    
    
    