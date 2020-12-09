function [] = Plot3D(Names,meas_type)

ct = length(Names);
CsV = hsv(ct);
h = figure();
set(h,'Position',[82 151 1497 827]);
for i = 1:ct-1
    [~,name,~] = fileparts(Names{i});
    PlotResults_3DSim_tmt(Names{i},CsV(i,:),name,h,meas_type);
end

end