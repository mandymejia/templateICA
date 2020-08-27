function [miu_svq_fixzq, miu_svq2_fixzq] = Cal_SvMmtCondzl( miu_svq_condz , miu_svq2_condz, PostProbz, idx, M, Q2, q2)
% [miu_svq_fixzq, miu_svq2_fixzq] = Cal_SvMmtCondzl( miu_svq_condz , miu_svq2_condz, PostProbz, idx, M, Q2, q2)
%
% INPUTS:
% miu_svq_condz (M^Q2 x 1)   : E[s_q(v)|z(v),Y(v),theta]
% miu_svq2_condz (M^Q2 x 1)  : E[s_q(v)^2|z(v),Y(v),theta]
% PostProbz (M^Q2 x 1) 		 : p(z(v)|Y(v),theta)
% idx (1 x M^(Q2-1) )  		 : indices of M^Q2 possible values of z(v) where z_q(v)=m
% M  (1x1) 			   		 : number of MoG components
% Q2 (1x1)			   		 : number of free ICs
% q2 (1x1)			   		 : index of free component (q2=1,...,Q2)
%
% OUTPUTS:
% miu_svq_fixzq          : first order moment of s_q(v) given z_q(v)=m (marginal)
% miu_svq2_fixzq         : second order moments of s_q(v) given z_q(v)=m (marginal)

miu_svq_fixzq  = 0;
miu_svq2_fixzq = 0;
ID = idx; 

sumprob = sum(PostProbz(idx)); %p(z_q(v)=m|Y(v),theta)
if(sumprob ~= 0)
    weight = PostProbz(idx)./sumprob ; % Prob(z_{-q}(v) | z_q(v) = m)
else
    weight = PostProbz(idx);
end;

%iteratively sum over all z(v) where z_q(v)=m ( M^(Q-1) possibilities)
for  ii = 1:M^(Q2-1)
     id = ID(ii); %indexes original M^Q2 possibilities
     miu_svq_fixzq  = miu_svq_fixzq  + miu_svq_condz(id)*weight(ii);
     miu_svq2_fixzq = miu_svq2_fixzq + miu_svq2_condz(id)*weight(ii);
end;
end










