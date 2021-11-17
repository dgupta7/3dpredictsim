function [] = PlotEnergySmoothVsNonSmooth(ResultsFile,LegName)



scs = get(0,'ScreenSize');
fsc = scs;
fsc(2) = 20;
fsc(4) = fsc(4)-120;

h = figure('Position',fsc);
h.Name = LegName;
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





load(ResultsFile,'R');

iM = [47:92];    

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
        
    subplot(5,4,ii)
    plot(R.MetabB.Etot(:,iM(i)),'DisplayName',['Smoothed (b=' num2str(R.S.tanh_b) ')']); hold on;
    plot(R.MetabB_non_smooth.Etot(:,iM(i)),'--','DisplayName','Non-smoothed')
    title(R.colheaders.muscles{iM(i)},'Interpreter','none')
    ylabel('Edot');
    if ii==4
        legend
    end
    
    subplot(5,4,ii+4)
    plot(R.MetabB.Adot(:,iM(i))); hold on;
    plot(R.MetabB_non_smooth.Adot(:,iM(i)),'--');
    ylabel('Adot');

    subplot(5,4,ii+8)
    plot(R.MetabB.Mdot(:,iM(i))); hold on;
    plot(R.MetabB_non_smooth.Mdot(:,iM(i)),'--');
    ylabel('Mdot');

    subplot(5,4,ii+12)
    plot(R.MetabB.Sdot(:,iM(i))); hold on;
    plot(R.MetabB_non_smooth.Sdot(:,iM(i)),'--');
    ylabel('Sdot');

    subplot(5,4,ii+16)
    plot(R.MetabB.Wdot(:,iM(i))); hold on;
    plot(R.MetabB_non_smooth.Wdot(:,iM(i)),'--');
    ylabel('Wdot');
    
    xlabel('% stride')

end




   
    
end
    
    
    