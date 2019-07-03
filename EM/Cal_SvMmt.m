function [miu_sv, miu_sv_svT] = Cal_SvMmt(miu_sv_condz , miu_sv_svT_condz, PostProbz, M, Q2)
% [miu_sv, miu_sv_svT] = Cal_SvMmt(miu_sv_condz , miu_sv_svT_condz, PostProbz, M)
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Input: 
% miu_sv_condz     - (Q x 1 x M^Q) E(s(v) | z(v), Y(v), theta) 
% miu_sv_svT_condz - (Q x Q x M^Q) E(s(v)s(v)' | z(v), Y(v), theta) 
% PostProbz        - (M^Q x 1)     Prob(z(v) | Y(v), theta) 
% where s(v) = [s1(v); s2(v); ...sQ(v)]    (Q x 1)
% Q2 - number of free components
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Output:
% E(s(v)      | Y(v), theta)
% E(s(v)s(v)' | Y(v), theta)

  Q = size(miu_sv_condz,1);
  miu_sv     = zeros(Q, 1); %first moments of posterior of s(v)
  miu_sv_svT = zeros(Q, Q); %second moments of posterior of s(v)
  %save all moments independent of z(v)
  for ii = 1:M^Q2

    %strcat(num2str(round(ii/(M^Q2)*100)),'%')

    if ii == 1
                  miu_sv     =  miu_sv_condz(:,:,ii)*PostProbz(ii);
                  miu_sv_svT =  miu_sv_svT_condz(:,:,ii)*PostProbz(ii);
    else
                  miu_sv      =  miu_sv +  miu_sv_condz(:,:,ii)*PostProbz(ii);
                  miu_sv_svT =  miu_sv_svT + miu_sv_svT_condz(:,:,ii)*PostProbz(ii);
    end;
  end;

end
