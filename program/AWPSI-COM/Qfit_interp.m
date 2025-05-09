function [S,S_med,S_med_cutnan_interp]=Qfit_interp(minQ,maxQ,window_size,S)

S(S > maxQ) = NaN;
S(S < minQ) = NaN;
S = inpaint_nans(S, 3); % 使用inpaint_nans函数进行插值
% figure;
% imagesc(S);
% colormap(jet);
% caxis([30 70]);

S_med=zeros(size(S));
for i = 1:size(S, 1)
    S_med(i, :) = medfilt1(S(i, :), window_size);
%         pi =polyfit(1:length(S(1, :)),S(i, :),30);
%     S_med(i, :) = polyval(pi,1:length(S(1, :)));
end
% figure;
% imagesc(S_med);
% colormap(jet);
% caxis([30 70]);

S_med_cutnan=S_med;
S_med_cutnan = S_med_cutnan(~any(isnan(S_med_cutnan), 2), :);
nrows = 2036;
orig_rows = linspace(1, size(S_med_cutnan, 1), size(S_med_cutnan, 1));
interp_rows = linspace(1, size(S_med_cutnan, 1), nrows);
S_med_cutnan_interp = zeros(nrows, size(S_med_cutnan, 2));
for j = 1:size(S_med_cutnan, 2)
    orig_col = S_med_cutnan(:, j);
    interp_col = interp1(orig_rows, orig_col, interp_rows, 'line');
    S_med_cutnan_interp(:, j) = interp_col;
end
% figure;
% imagesc(B);
% colormap(jet);
% caxis([30 70]);
% 
% figure;
% imagesc(Q_eff);
% colormap('jet');
% colorbar;
% caxis([30 70]);