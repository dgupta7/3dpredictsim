function [MA_PF, l_PF, v_PF] = getPlantarFasciaLengthVelocity(q_mtj,qdot_mtj,q_mtp,qdot_mtp,varargin)

if isempty(varargin)
    R_mtth = 9.5e-3; % radius of the first metatarsal head
else
    R_mtth = varargin{1};
end

%%
% see \FootModel\windlassGeometryPolynomials.m
% length of PF spanning arch
l_PF_fa = 0.1392179 + 0.0374482.*q_mtj.^1 + -0.0166876.*q_mtj.^2 + ...
    -0.001758651.*q_mtj.^3 + 0.0004480769.*q_mtj.^4;
% moment arm of PF to mtj
MA_PF = 0.0374478 + -0.03337403.*q_mtj.^1 + -0.005255987.*q_mtj.^2 + ...
        0.001767266.*q_mtj.^3 + -0.0001071423.*q_mtj.^4 + 9.858065e-05.*q_mtj.^5;

% Plantar fascia length
l_PF = l_PF_fa + R_mtth*(pi/4+q_mtp) + 0.0042; % constant term to correct to physical length

%%

v_PF = (0.0374482 + -0.0166876.*q_mtj*2 + -0.001758651.*q_mtj.^2*3 +...
    0.0004480769.*q_mtj.^3*4).*qdot_mtj + R_mtth*qdot_mtp;




end