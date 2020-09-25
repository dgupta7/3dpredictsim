%% Test time shift exoskeleton assistance

pathRepo = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim';
Poggensee = load([pathRepo,'\Data\Poggensee_2020_Exp\torque_profile.mat']);
Tankle = Poggensee.torque; % -1 because plantarflexion is negative in opensim model
nfr = length(Tankle);

figure();
for xNew = 0.5:0.01:0.7
        
    Tankle_stance = Tankle(1:nfr*0.61);
    Tankle_swing = Tankle(nfr*0.61:nfr);
    
    SplineStance = spline(linspace(0,1,length(Tankle_stance)),Tankle_stance');
    SplineSwing = spline(linspace(0,1,length(Tankle_swing)),Tankle_swing');
    xStance = linspace(0,1,nfr*xNew);
    xSwing = linspace(0,1,nfr*(1-xNew));
    
    TStance = ppval(SplineStance,xStance);
    TSwing = ppval(SplineSwing,xSwing);
    Tnew = [TStance TSwing(2:end)];    
    
    plot(Tnew,'r');  hold on;
end
plot(Tankle,'b'); hold on;



