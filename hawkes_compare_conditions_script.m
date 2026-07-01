cd     '/Users/wehr/Documents/Analysis/geotrig'



c1=load('pp_analysis-output-29-Jun-2026-18:48:15.mat')
c1.condition_name
%dark laseron

c2=load('pp_analysis-output-29-Jun-2026-18:57:10.mat');
c2.condition_name
%light laseron

c3=load('pp_analysis-output-29-Jun-2026-17:34:41.mat');
%pp_hawkes_plot( eventnames, c2.hawkes_params, c2.condition_name)
c3.condition_name
%dark laseroff

c4=load('pp_analysis-output-29-Jun-2026-17:23:50.mat')
c4.condition_name
%pp_hawkes_plot( eventnames, c1.hawkes_params, c1.condition_name)
%light laseroff


figure
diffMat=c2.hawkes_params.branchingRatio' - c1.hawkes_params.branchingRatio';
mytitle=sprintf('Difference hawkes, %s - %s', c2.condition_name, c1.condition_name );
plot_tpm_diff_circle(diffMat, eventnames, 'Title', mytitle)

C={c1, c2, c3, c4};
figure;
tiledlayout(2,2, "TileSpacing","compact")
set(gcf, "Position", [460 460 1020 840])
% 1-3 : dark on-off
% 2-4 : light on-off
% 1-2 : dark-light on
% 3-4 : dark-light off
A=[1 2 1 3];
B=[3 4 2 4];

for c=1:4
    a=A(c);
    b=B(c);
    nexttile

    diffMat=C{a}.hawkes_params.branchingRatio' - C{b}.hawkes_params.branchingRatio';
    mytitle=sprintf('Difference hawkes, %s - %s', C{a}.condition_name, C{b}.condition_name );
    plot_tpm_diff_circle(diffMat, eventnames, 'Title', mytitle)

end


