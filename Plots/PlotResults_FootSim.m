function [varargout] = PlotResults_FootSim(R,varargin)
fig2 = 0;

% if ~isfield(R,'legname')
%     if ~isfield(R.S,'mtj_stiffness')
%         R.S.mtj_stiffness = '';
%     end
%     legname = ['PF: ' R.PF_stiffness ', l_s = ' num2str(R.S.PF_slack_length*1000) '; mtj: ' R.S.mtj_stiffness];
%     R.legname = legname;
% end

n_mtp = length(R.Qs_mtp);
n_tib = length(R.Fs_tib);

mrk = {'-o','--d','-.v',':*','-p','-x','-s','-h','-+','-<','-^'};

BW = R.S.mass*9.81;

if ~isempty(varargin)
    CsV = varargin{1};
else
    CsV = 'r';
end
if length(varargin) >= 2
    h = varargin{2};
    tab1 = h.Parent.Children(1).Children(1).Children(1);
    tab2 = h.Parent.Children(1).Children(1).Children(2);
    tab3 = h.Parent.Children(1).Children(1).Children(3);
    tab4 = h.Parent.Children(1).Children(1).Children(4);
    tab5 = h.Parent.Children(1).Children(1).Children(5);
    tab6 = h.Parent.Children(1).Children(1).Children(6);
    BoolFirst = 0;
else
    h = figure('Position',[82,151,1497,827]);
    h.Name = 'FootModel_Results';
    hTabGroup = uitabgroup;
    tab1 = uitab(hTabGroup, 'Title', 'Unloaded');
    tab2 = uitab(hTabGroup, 'Title', 'Load effects 1');
    tab3 = uitab(hTabGroup, 'Title', 'Load effects 2');
    tab4 = uitab(hTabGroup, 'Title', 'Arch stiffness');
    tab5 = uitab(hTabGroup, 'Title', 'mtp');
    tab6 = uitab(hTabGroup, 'Title', 'mtp stiffening');
    BoolFirst = 1;
end

% return figure handle if wanted
if nargout == 1
    varargout{1} = h;
end

%% unloaded condition

axes('parent', tab1);


subplot(3,3,1)
hold on
if length(R.Qs_mtp)>1 && min(R.Qs_mtp)*180/pi==-20
    l0 = max(R.l_fa(:,1));
    plot(R.Qs_mtp*180/pi,R.l_fa(:,1)./l0,'.','color',CsV)
    c = polyfit(R.Qs_mtp'*180/pi,R.l_fa(:,1)./l0,1);
    lm = polyval(c,R.Qs_mtp*180/pi);
    plot(R.Qs_mtp*180/pi,lm,'-','color',CsV)
    
    ylabel('arch length normalised to max (-)')
    
else
    plot(R.Qs_mtp*180/pi,R.l_fa(:,1)*1000,'color',CsV)
    ylabel('arch length (mm)')
    
end
xlabel('mtp angle (°)')
title('Foot arch length')


subplot(3,3,2)
hold on
plot(R.Qs_mtp*180/pi,R.h_fa(:,1)*1000,'color',CsV)
xlabel('mtp angle (°)')
ylabel('arch height (mm)')
title('Foot arch height')

subplot(3,3,3)
hold on
plot(R.Qs_mtp*180/pi,R.Qs(:,1,R.jointfi.tmt.r)*180/pi,'color',CsV,'DisplayName',R.legname)
xlabel('mtp angle (°)')
ylabel('mt angle (°)')
title('Midtarsal joint angle')
lg12 = legend('Location','northeast');
lhPos = lg12.Position;
lhPos(1) = lhPos(1)+0.1;
set(lg12,'position',lhPos);

subplot(3,3,4)
hold on
plot(R.Qs_mtp*180/pi,R.l_PF(:,1)*1000,'color',CsV)
xlabel('mtp angle (°)')
ylabel('PF length (mm)')
title('Plantar fascia length')

subplot(3,3,5)
hold on
if isfield(R,'MA_PF')
    plot(R.Qs_mtp*180/pi,R.MA_PF(:,1)*1000,'color',CsV)
end
xlabel('mtp angle (°)')
ylabel('PF moment arm (mm)')
title('Plantar fascia moment arm')

subplot(3,3,6)
hold on
plot(R.Qs_mtp*180/pi,R.M_li(:,1),'color',CsV)
xlabel('mtp angle (°)')
ylabel('Torque (Nm)')
title('mt torque (except PF)')

subplot(3,3,7)
hold on
plot(R.Qs_mtp*180/pi,R.F_PF(:,1),'color',CsV)
xlabel('mtp angle (°)')
ylabel('F PF (N)')
title('Plantar fascia force')

subplot(3,3,8)
hold on
plot(R.Qs_mtp*180/pi,R.M_PF(:,1),'color',CsV)
xlabel('mtp angle (°)')
ylabel('Torque (Nm)')
title('Plantar fascia torque')

subplot(3,3,9)
hold on
plot(R.Qs_mtp*180/pi,R.GRF_calcn(:,1,2),'-','color',CsV)
hold on
plot(R.Qs_mtp*180/pi,R.GRF_metatarsi(:,1,2),'-.','color',CsV)
xlabel('mtp angle (°)')
ylabel('GRF_y (N)')
title('vertical GRF')
lg13=legend('calcaneus','metatarsi','Location','northeast');
lhPos = lg13.Position;
lhPos(1) = lhPos(1)+0.1;
set(lg13,'position',lhPos);

%%

for i=1:n_mtp
    axes('parent', tab2);
    
    js = find(R.failed(i,:)==0);
    Fs_tib = R.Fs_tib(js);
    
    subplot(2,3,1)
    hold on
    plot(Fs_tib,R.l_fa(i,js)*1000,mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('arch length (mm)')
    title('Foot arch length (sagittal)')

    subplot(2,3,2)
    hold on
    plot(Fs_tib,R.h_fa(i,js)*1000,mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('arch height (mm)')
    title('Foot arch height (sagittal)')
    
    subplot(2,3,3)
    hold on
    if isfield(R,'MA_PF')
        plot(Fs_tib,R.MA_PF(i,js)*1000,mrk{i},'Color',CsV)
    end
    xlabel('tibia force (N)')
    ylabel('PF moment arm (mm)')
    title('Plantar fascia moment arm')
    
    subplot(2,3,4)
    hold on
    plot(Fs_tib,R.l_PF(i,js)*1000,mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('PF length (mm)')
    title('Plantar fascia length (sagittal)')
    
    subplot(2,3,5)
    hold on
    plot(Fs_tib,R.F_PF(i,js),mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Qs_mtp(i)*180/pi) '°; ' R.legname])
    xlabel('tibia force (N)')
    ylabel('F PF (N)')
    title('Plantar fascia force')
    lg02=legend('Location','northeast');

end

title(lg02,'mtp, model')
lhPos = lg02.Position;
lhPos(1) = lhPos(1)+0.25;
set(lg02,'position',lhPos);

%%
for i=1:n_mtp
    axes('parent', tab3);
    
    js = find(R.failed(i,:)==0);
    Fs_tib = R.Fs_tib(js);

    subplot(3,3,1)
    hold on
    plot(Fs_tib,R.Qs(i,js,R.jointfi.ankle.r)*180/pi,mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('angle (°)')
    title('ankle angle')

    subplot(3,3,2)
    hold on
    plot(Fs_tib,R.Qs(i,js,R.jointfi.subt.r)*180/pi,mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('angle (°)')
    title('subtalar angle')
    
    subplot(3,3,3)
    hold on
    plot(Fs_tib,R.Qs(i,js,R.jointfi.tmt.r)*180/pi,mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Qs_mtp(i)*180/pi) '°, ' R.legname])
    xlabel('tibia force (N)')
    ylabel('angle (°)')
    lg02=legend('Location','northeast');
    title('midtarsal angle')
    
    subplot(3,3,4)
    hold on
    plot(Fs_tib,R.T_ankle.ext(i,js),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('torque (Nm)')
    title('ankle torque')

    subplot(3,3,5)
    hold on
    plot(Fs_tib,R.T_subt.ext(i,js),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('torque (Nm)')
    title('subtalar torque')
    
    subplot(3,3,6)
    hold on
    plot(Fs_tib,R.M_li(i,js),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('torque (Nm)')
    title('midtarsal torque (except PF)')

    subplot(3,3,7)
    hold on
    plot(Fs_tib,R.GRF_calcn(i,js,2),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('GRF_y (N)')
    title('vertical GRF calcaneus')

    subplot(3,3,8)
    hold on
    plot(Fs_tib,R.GRF_metatarsi(i,js,2),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('GRF_y (N)')
    title('vertical GRF metatarsi')
    
    subplot(3,3,9)
    hold on
    plot(Fs_tib,R.M_PF(i,js),mrk{i},'Color',CsV)
    xlabel('tibia force (N)')
    ylabel('torque (Nm)')
    title('midtarsal torque (PF)')
    
end

title(lg02,'mtp, model')
lhPos = lg02.Position;
lhPos(1) = lhPos(1)+0.1;
set(lg02,'position',lhPos);

%%
axes('parent', tab5);

for i=1:n_tib
    
    js = find(R.failed(:,i)==0);
    
    subplot(2,3,1)
    hold on
    plot(R.Qs_mtp(js)*180/pi,R.Qs(js,i,R.jointfi.mtp.r)*180/pi,mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Fs_tib(i)) 'N, ' R.legname])
    xlabel('angle wrt ground(°)')
    ylabel('angle wrt forefoot(°)')
    title('mtp angle')
    
    subplot(2,3,2)
    hold on
    plot(R.Qs_mtp(js)*180/pi,R.Qs(js,i,R.jointfi.tmt.r)*180/pi,mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Fs_tib(i)) 'N, ' R.legname])
    xlabel('mtp angle wrt ground(°)')
    ylabel('angle(°)')
    title('mtj angle')
    
    subplot(2,3,3)
    hold on
    T_mtp = R.M_mtp(js,i);% - R.S.kMTP*R.Qs_mtp(:);
    plot(R.Qs_mtp(js)*180/pi,T_mtp,mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Fs_tib(i)) 'N, ' R.legname])
    xlabel('mtp angle wrt ground(°)')
    ylabel('torque(Nm)')
    title('mtp torque from PF')
    lg03=legend('Location','northeast');
    
    subplot(2,3,4)
    hold on
    plot(R.Qs_mtp(js)*180/pi,R.l_PF(js,i)*1000,mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Fs_tib(i)) 'N, ' R.legname])
    xlabel('mtp angle wrt ground(°)')
    ylabel('length (mm)')
    title('PF length')
    
    if isfield(R,'T_mtp')
        subplot(2,3,5)
        hold on
        plot(R.Qs_mtp(js)*180/pi,R.T_mtp(js,i),mrk{i},'Color',CsV,...
        'DisplayName',[num2str(R.Fs_tib(js(i))) 'N, ' R.legname])
        xlabel('mtp angle wrt ground(°)')
        ylabel('torque(Nm)')
        title('external mtp torque')
    end
end
title(lg03,'F, model')
lhPos = lg03.Position;
lhPos(2) = lhPos(2)-0.4;
set(lg03,'position',lhPos); 

    
%%
axes('parent', tab6);

for i=1:n_tib
    pl = polyfit(R.Qs_mtp(js)'*180/pi,T_mtp,1);
    pln = polyfit(R.Qs(js,i,R.jointfi.mtp.r)*180/pi,T_mtp,1);
    
    kg(i) = -pl(1);
    kf(i) = -pln(1);

    if isfield(R,'T_mtp')
        ple = polyfit(R.Qs_mtp(js)'*180/pi,R.T_mtp(js,i),1);
        plne = polyfit(R.Qs(js,i,R.jointfi.mtp.r)*180/pi,R.T_mtp(js,i),1);
        
        kge(i) = ple(1);
        kfe(i) = plne(1);
    end

end

subplot(2,4,1)
hold on
grid on
% p11=plot(Fs_tib(js)/BW,kg,'o','Color',CsV,'DisplayName','toes wrt ground');
p12=plot(R.Fs_tib/BW,kf,'d','Color',CsV,'DisplayName','toes wrt forefoot');
xlabel('vertical tibia load / BW (-)')
ylabel('quasi-stiffness (Nm/°)')
title('quasi-stiffness mtp (from PF only)')
% lg05=legend([p11,p12]);
% title(lg05,'mtp angle definition')
% lhPos = lg05.Position;
% lhPos(2) = lhPos(2)-0.1;
% set(lg05,'position',lhPos); 

subplot(2,4,2)
hold on
grid on
% plot(Fs_tib(js)/BW,kg./kg(1),'o','Color',CsV,'DisplayName',R.legname)
plot(R.Fs_tib/BW,kf./kf(1),'d','Color',CsV,'DisplayName',R.legname)
xlabel('vertical tibia load / BW (-)')
ylabel('relative quasi-stiffness (-)')
title('quasi-stiffness mtp (from PF only)')
lg06=legend('location','northwest');
title(lg06,'plantar fascia and midtarsal joint stiffness models')
lhPos = lg06.Position;
lhPos(2) = lhPos(2)-0.4;
set(lg06,'position',lhPos); 


if isfield(R,'T_mtp')
    subplot(2,4,3)
    hold on
    grid on
%     plot(Fs_tib(js)/BW,kge,'o','Color',CsV,'DisplayName',R.legname)
    plot(R.Fs_tib/BW,kfe,'d','Color',CsV,'DisplayName',R.legname)
    xlabel('tibia force / BW (-)')
    ylabel('quasi-stiffness (Nm/°)')
    title('quasi-stiffness mtp (external)')

    subplot(2,4,4)
    hold on
    grid on
%     plot(Fs_tib(js)/BW,kge./kge(1),'o','Color',CsV,'DisplayName',R.legname)
    plot(R.Fs_tib/BW,kfe./kfe(1),'d','Color',CsV,'DisplayName',R.legname)
    xlabel('tibia force / BW (-)')
    ylabel('relative quasi-stiffness (-)')
    title('quasi-stiffness mtp (external)')
    
    subplot(2,4,7)
    hold on
    grid on
    pln = polyfit(R.Fs_tib/BW,kfe./kfe(1),1);
    plot(R.S.kMT_li,pln(1),'d','Color',CsV,'DisplayName',R.legname)
    xlabel('mtj stiffness (Nm/rad)')
    ylabel('mtp quasi-stiffening with load')
    title('mtp quasi-stiffness increase with load')
end

%%
axes('parent', tab4);

j = find(R.Qs_mtp(:)==0);

js = find(R.failed(j,:)==0);
Fs_tib = R.Fs_tib(js);
l_fa = R.l_fa(j,js);
h_fa = R.h_fa(j,js);

if BoolFirst
    % load reference graph images
    pathmain        = pwd;
    [pathRepo,~,~]  = fileparts(pathmain);
    folder = '\Figures';
    file = 'arch_stiffness_Ker87.png';
    pathRefImg = fullfile(pathRepo,folder,file);
    img_Ker = imread(pathRefImg);
    file = 'arch_stiffness_Welte18.png';
    pathRefImg = fullfile(pathRepo,folder,file);
    img_Welte = imread(pathRefImg);
    file = 'arch_stiffness_Stearne16.png';
    pathRefImg = fullfile(pathRepo,folder,file);
    img_Stearne = imread(pathRefImg);

end

subplot(3,5,[1,2,6,7])
if R.F_PF(j,js(end))<1
    plot((l_fa-R.L0)*1000,Fs_tib/1000,'color',CsV,'DisplayName',R.legname)
else
    plot((l_fa-l_fa(1))*1000,Fs_tib/1000,'color',CsV,'DisplayName',R.legname)
end
hold on
if BoolFirst
    hi1 = image([1,9.65]*0.9,flip([0,4]*0.9^2),img_Ker);
    uistack(hi1,'bottom')
end
xlabel('horizontal elongation (mm)')
ylabel('vertical force (kN)')
title({'Foot arch stiffness','as defined by Ker et al, 1987'})
lg12 = legend('Location','southeast');
lhPos = lg12.Position;
lhPos(2) = lhPos(2)-0.3;
set(lg12,'position',lhPos);



subplot(3,5,[5,10])
if R.F_PF(j,js(end))<1
    plot((R.H0-h_fa)*1000,Fs_tib/1000,'color',CsV)
else
    plot((h_fa(1)-h_fa)*1000,Fs_tib/1000,'color',CsV)
end
hold on
axis tight
if BoolFirst
    hi2 = image([-2,12],flip([-0.3,6]),img_Stearne);
    uistack(hi2,'bottom')
end
xlabel('vertical displacement (mm)')
ylabel('vertical force (kN)')
title({'Foot arch stiffness','as defined by Stearne et al, 2016'})

subplot(3,5,15)
if R.F_PF(j,js(end))<1
    plot((R.H0-h_fa(Fs_tib<=300))*1000,Fs_tib(Fs_tib<=300),'color',CsV)
else
    plot((h_fa(1)-h_fa(Fs_tib<=300))*1000,Fs_tib(Fs_tib<=300),'color',CsV)
end
hold on
xlabel('vertical displacement (mm)')
ylabel('vertical force (N)')
title('Detail view')


j1 = find(R.Qs_mtp(:)==-30*pi/180);
j2 = find(R.Qs_mtp(:)==30*pi/180);

if ~isempty(j1) && ~isempty(j2)
    idx = [j1,j2];
    n_i = 2;
else
    idx = 1:n_mtp;
    n_i = n_mtp;
end

for i=1:n_i
    
    idx_ac = find(R.failed(idx(i),:)==0 & R.Fs_tib<=BW);
    F_ac{i} = R.Fs_tib(idx_ac);
    h0_ac = R.h_fa(idx(i),1);
    tmp = h0_ac - R.h_fa(idx(i),idx_ac);
    ac{i} = tmp;
    ac_max(i) = max(ac{i});
end

subplot(3,5,[3,4,8,9])
hold on
for i=1:n_i
    ac_rel = ac{i}/max(ac_max);
    plot(ac_rel(:),F_ac{i}/BW,mrk{i},'color',CsV,'DisplayName',['mtp: ' num2str(R.Qs_mtp(idx(i))*180/pi) '°'])
end
if BoolFirst
    hi3 = image([0,1],flip([0,0.9]),img_Welte);
    uistack(hi3,'bottom')
end
xlabel('arch compression (-)')
ylabel('vertical force / body weight (-)')
title({'Foot arch stiffness','as defined by Welte et al, 2018'})
leg = legend('Location','southeast');
title(leg,'mtp dorsiflexion')
xlim([0,1])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if fig2
    % 
    h2 = figure('Position',[82,151,1497,827]);
    h2.Name = R.PF_stiffness;
    CsV1 = hsv(n_mtp);
    CsV2 = hsv(n_tib);

    %%

    for i=1:n_mtp

        figure(h2)
        subplot(2,3,[1,4])
        xy = [squeeze(R.toes_or(i,1,1:2)),squeeze(R.metatarsi_or(i,1,1:2)),...
            squeeze(R.calcn_or(i,1,1:2)),squeeze(R.talus_or(i,1,1:2)),...
            squeeze(R.tibia_or(i,1,1:2))];

        hold on
        grid on
        plot(xy(1,:),xy(2,:),'-o','Color',CsV1(i,:),'DisplayName',['mtp: ' num2str(R.Qs_mtp(i)*180/pi) '°'])
        axis equal

    end

    title('sagittal plane (right)')
    xlabel('x')
    ylabel('y')
    lg_pos=legend('Location','best');
    title(lg_pos,'F_y = 0N')

    %%

    j = find(R.Qs_mtp(:)==0);
    leg1 = [];

    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
    zmin = 0;
    zmax = 0;

    for i=1:(ceil(n_tib/5)):n_tib
        if R.failed(j,i) == 0
            figure(h2)

            xy = [squeeze(R.toes_or(j,i,1:2)),squeeze(R.metatarsi_or(j,i,1:2)),...
                squeeze(R.calcn_or(j,i,1:2)),squeeze(R.talus_or(j,i,1:2)),...
                squeeze(R.tibia_or(j,i,1:2))];

            z = [squeeze(R.toes_or(j,i,3)),squeeze(R.metatarsi_or(j,i,3)),...
                squeeze(R.calcn_or(j,i,3)),squeeze(R.talus_or(j,i,3)),...
                squeeze(R.tibia_or(j,i,3))];

            xmin = min(xmin,min(xy(1,:)));
            xmax = max(xmax,max(xy(1,:)));
            ymin = min(ymin,min(xy(2,:)));
            ymax = max(ymax,max(xy(2,1:end-1)));
            zmin = min(zmin,min(z));
            zmax = max(zmax,max(z));

            range_x = [xmin-0.01,xmax+0.01];
            range_y = [ymin-0.01,ymin-0.01+norm(range_x)];

            ankle_axis = squeeze(R.talus_or(j,i,:)) + [-1,1].*squeeze(R.ankle_axis(j,i,:))/100;
            COP_calcn = squeeze(R.COP_calcn(j,i,[1,3]));
            COP_calcn(2) = -COP_calcn(2); %to plot -z
            COP_metatarsi = squeeze(R.COP_metatarsi(j,i,[1,3]));
            COP_metatarsi(2) = -COP_metatarsi(2);
            COP_R = [R.GRF_calcn(j,i,2),R.GRF_metatarsi(j,i,2)];
            COP_R = COP_R./sum(COP_R)/30;
            
            subplot(2,3,2)
            hold on
            plot(xy(1,4:5),xy(2,4:5),'-','Color',CsV2(i,:))
            plot(xy(1,1),xy(2,1),'x','Color',CsV2(i,:))
            plot(xy(1,2),xy(2,2),'d','Color',CsV2(i,:))
            plot(xy(1,3),xy(2,3),'.','Color',CsV2(i,:))
            plot(xy(1,4),xy(2,4),'*','Color',CsV2(i,:))
            plot(ankle_axis(1,:),ankle_axis(2,:),'--','Color',CsV2(i,:))
            axis equal
            title('sagittal plane (right)')
            xlabel('x')
            ylabel('y')
            if i>1
                ylim(range_y)
                xlim(range_x)
            end
            range_y = get(gca, 'ylim');

            subplot(2,3,3)
            hold on
            plot(-z(1,4:5),xy(2,4:5),'-','Color',CsV2(i,:))
            p1=plot(-z(1,1),xy(2,1),'x','Color',CsV2(i,:),'DisplayName',['F_y = ' num2str(R.Fs_tib(i)) 'N']);
            plot(-z(1,2),xy(2,2),'d','Color',CsV2(i,:))
            plot(-z(1,3),xy(2,3),'.','Color',CsV2(i,:))
            plot(-z(1,4),xy(2,4),'*','Color',CsV2(i,:))
            plot(-ankle_axis(3,:),ankle_axis(2,:),'--','Color',CsV2(i,:))
            leg1(end+1) = p1;
            axis equal
            title('frontal plane (front)')
            xlabel('-z')
            ylabel('y')
            if i>1
                ylim(range_y)
            end
            range_z = get(gca, 'xlim');
            lg1=legend(leg1,'Location','northeast');

            subplot(2,3,5)
            hold on
            p6=plot(xy(1,4:5),-z(1,4:5),'-','Color',CsV2(i,:),'DisplayName','tibia');
            p2=plot(xy(1,1),-z(1,1),'x','Color',CsV2(i,:),'DisplayName','toes or');
            p3=plot(xy(1,2),-z(1,2),'d','Color',CsV2(i,:),'DisplayName','metatarsi or');
            p4=plot(xy(1,3),-z(1,3),'.','Color',CsV2(i,:),'DisplayName','calcn or');
            p5=plot(xy(1,4),-z(1,4),'*','Color',CsV2(i,:),'DisplayName','talus or');
            p7=plot(ankle_axis(1,:),-ankle_axis(3,:),'--','Color',CsV2(i,:),'DisplayName','ankle axis');
%             p8=plot(COP(1,2:3),-COP(3,2:3),'o','Color',CsV2(i,:),'DisplayName','centres of pressure');
            if i>1
                p8=viscircles(COP_calcn',COP_R(1),'Color',CsV2(i,:),'LineWidth',1,'LineStyle',':');
                p9=viscircles(COP_metatarsi',COP_R(2),'Color',CsV2(i,:),'LineWidth',1,'LineStyle',':');
                set(p8,'DisplayName','COP calcn')
                set(p9,'DisplayName','COP metatarsi')
                lg2=legend([p2,p3,p4,p5,p6,p7,p8,p9],'Location','northeast');
            else
                lg2=legend([p2,p3,p4,p5,p6,p7],'Location','northeast');
            end
            axis equal
            title('transverse plane (top)')
            xlabel('x')
            ylabel('-z')
            if i>1
                xlim(range_x)
                ylim(range_z)
            end
            
        end
    end


    figure(h2)
    title(lg1,'colour meaning')
    lhPos = lg1.Position;
    lhPos(2) = lhPos(2)-0.4;
    % lhPos(1) = lhPos(1)+0.1;
    set(lg1,'position',lhPos);

    title(lg2,'symbol meaning')
    lhPos = lg2.Position;
    lhPos(1) = lhPos(1)+0.2;
    set(lg2,'position',lhPos);




    figure(h)
end

end