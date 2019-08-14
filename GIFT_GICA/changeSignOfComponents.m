function [A, W, icasig, skew] = changeSignOfComponents(A, icasig)
% Change sign of components

numOfIC = size(icasig, 1);
skew = zeros(1, numOfIC);
%force group images to be positive
for compNum = 1:numOfIC
    v = icatb_recenter_image(icasig(compNum, :));
    skew(compNum) = icatb_skewness(v) + eps;
    clear v;
    if (sign(skew(compNum)) == -1)
        %disp(['Changing sign of component ',num2str(compNum)]);
        icasig(compNum, :) = icasig(compNum, :)*-1;
        A(:, compNum) = A(:, compNum)*-1;
    end
    
end

W = pinv(A);
