function [] = Plot3D(Names,meas_type,varargin)

if length(varargin)>=1
    if varargin{1} == 0
        plot_normal = 1;
        plot_val = 0;
    elseif varargin{1} == 1
        plot_normal = 0;
        plot_val = 1;
    elseif varargin{1} == 2
        plot_normal = 1;
        plot_val = 1;
    else
        plot_normal = 0;
        plot_val = 0;
    end
else
    plot_normal = 1;
    plot_val = 0;
end

ct = length(Names);
CsV = hsv(ct);
if plot_normal
    h = figure();
    set(h,'Position',[82 151 1497 827]);
    for i = 1:ct
        [~,name,~] = fileparts(Names{i});
        PlotResults_3DSim_tmt(Names{i},CsV(i,:),name,h,meas_type);
    end
end

if plot_val
    h1 = figure();
    set(h1,'Position',[82 151 1497 827]);
    for i = 1:ct
        [~,name,~] = fileparts(Names{i});
        ValidationPlots(Names{i},CsV(i,:),name,h1);
    end
end

end