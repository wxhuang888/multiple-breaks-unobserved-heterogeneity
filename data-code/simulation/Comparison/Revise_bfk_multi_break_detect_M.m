function [est_break_pts] = Revise_bfk_multi_break_detect_M(N,T,ymat,xmat,zmat,ifull, M)

% This function is adpated from bfk_multi_break_detect.m
% The # of break M is an input and known (M = 1,2, or 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimate three structural breaks using sequential least 
% square on heterogeneous panels proposed by Baltagi, Feng, and Kao (2016)
%
% Input: 
% N: cross-sectional dimension
% T: time dimension
% ymat: y in a matrix form 
% xmat: x in a matrix form
% zmat: x whose coefficients contain breaks
% ifull: indicator whether the whole coefficient vector contain breaks
% 
% Output: 
% est_break_pts: estimated break points
% bhat: estimated coefficients in each regime
% sehat: estimated standard error correpsonding to bhat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if M > 3 | M < 1;
    error('check M!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate the first break    
K = size(xmat,3);
[est_break_pts(1)] = bfk_single_break_detect(N,T,xmat,ymat,zmat);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if M > 1
%% estimate the second break
if est_break_pts(1) <= 2 
    y_pd2 = ymat(est_break_pts(1)+1:T,:);
    x_pd2 = xmat(est_break_pts(1)+1:T,:,:);
    if ifull == 1
        z_pd2 = x_pd2(:,:,1:end);
    else
        z_pd2 = x_pd2(:,:,end);
    end
    [breakpt_tmp2] = bfk_single_break_detect(N,T-est_break_pts(1),x_pd2,y_pd2,z_pd2);
    est_break_pts(2) = est_break_pts(1)+breakpt_tmp2;
elseif est_break_pts(1) >= T-2  
    y_pd1 = ymat(1:est_break_pts(1),:);
    x_pd1 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd1 = x_pd1(:,:,1:end);
    else
        z_pd1 = x_pd1(:,:,end);
    end    
    [breakpt_tmp1] = bfk_single_break_detect(N,est_break_pts(1),x_pd1,y_pd1,z_pd1);
    est_break_pts(2) = breakpt_tmp2;
else    
    y_pd1 = ymat(1:est_break_pts(1),:);
    x_pd1 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd1 = x_pd1(:,:,1:end);
    else
        z_pd1 = x_pd1(:,:,end);    
    end      
    [breakpt_tmp1] = bfk_single_break_detect(N,est_break_pts(1),x_pd1,y_pd1,z_pd1);

    y_pd2 = ymat(est_break_pts(1)+1:T,:);
    x_pd2 = xmat(est_break_pts(1)+1:T,:,:);
    if ifull == 1
        z_pd2 = x_pd2(:,:,1:end);
    else
        z_pd2 = x_pd2(:,:,end);
    end
    [breakpt_tmp2] = bfk_single_break_detect(N,T-est_break_pts(1),x_pd2,y_pd2,z_pd2);
    
    y_pd11_tmp = ymat(1:breakpt_tmp1,:);
    x_pd11_tmp = xmat(1:breakpt_tmp1,:,:);
    [ssr11,bhat11] = bfk_post_break_est(N,breakpt_tmp1,y_pd11_tmp,x_pd11_tmp);   
    y_pd12_tmp = ymat(breakpt_tmp1+1:est_break_pts(1),:);
    x_pd12_tmp = xmat(breakpt_tmp1+1:est_break_pts(1),:,:);
    [ssr12,bhat12] = bfk_post_break_est(N,est_break_pts(1)-breakpt_tmp1,y_pd12_tmp,x_pd12_tmp);
    y_pd13_tmp = ymat(est_break_pts(1)+1:T,:);
    x_pd13_tmp = xmat(est_break_pts(1)+1:T,:,:);
    [ssr13,bhat13] = bfk_post_break_est(N,T-est_break_pts(1),y_pd13_tmp,x_pd13_tmp);
        
    y_pd21_tmp = ymat(1:est_break_pts(1),:);
    x_pd21_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr21,bhat21] = bfk_post_break_est(N,est_break_pts(1),y_pd21_tmp,x_pd21_tmp);   
    y_pd22_tmp = ymat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp2,:);
    x_pd22_tmp = xmat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp2,:,:);
    [ssr22,bhat22] = bfk_post_break_est(N,breakpt_tmp2,y_pd22_tmp,x_pd22_tmp);
    y_pd23_tmp = ymat(est_break_pts(1)+breakpt_tmp2+1:T,:);
    x_pd23_tmp = xmat(est_break_pts(1)+breakpt_tmp2+1:T,:,:);
    [ssr23,bhat23] = bfk_post_break_est(N,T-est_break_pts(1)-breakpt_tmp2,y_pd23_tmp,x_pd23_tmp);

    ssr1 = ssr11 + ssr12 + ssr13;
    ssr2 = ssr21 + ssr22 + ssr23;    
    if ssr1 < ssr2
        est_break_pts(2) = breakpt_tmp1;
    else
        est_break_pts(2) = est_break_pts(1)+breakpt_tmp2;
    end
end
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if M > 2
%% estimate the third break
est_break_pts = sort(est_break_pts);
if est_break_pts(1) <= 2 && est_break_pts(2)-est_break_pts(1)<=2
    y_pd3 = ymat(est_break_pts(2)+1:T,:);
    x_pd3 = xmat(est_break_pts(2)+1:T,:,:);
    if ifull == 1
        z_pd3 = x_pd3(:,:,1:end);        
    else
        z_pd3 = x_pd3(:,:,end);        
    end    
    [breakpt_tmp3] = bfk_single_break_detect(N,T-est_break_pts(2),x_pd3,y_pd3,z_pd3);
    est_break_pts(3) = est_break_pts(2)+breakpt_tmp3;
    
elseif est_break_pts(2) >= T-2 && est_break_pts(2)-est_break_pts(1)<=2
    y_pd3 = ymat(1:est_break_pts(1),:);
    x_pd3 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd3 = x_pd3(:,:,1:end);        
    else
        z_pd3 = x_pd3(:,:,end);        
    end
    [breakpt_tmp3] = bfk_single_break_detect(N,est_break_pts(1),x_pd3,y_pd3,z_pd3);
    est_break_pts(3) = breakpt_tmp3;

elseif est_break_pts(1) <= 2
    y_pd31 = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd31 = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    if ifull == 1
        z_pd31 = x_pd31(:,:,1:end);
    else
        z_pd31 = x_pd31(:,:,end);
    end       
    [breakpt_tmp31] = bfk_single_break_detect(N,est_break_pts(2)-est_break_pts(1),x_pd31,y_pd31,z_pd31);

    y_pd32 = ymat(est_break_pts(2)+1:T,:);
    x_pd32 = xmat(est_break_pts(2)+1:T,:,:);
    if ifull == 1
        z_pd32 = x_pd32(:,:,1:end);
    else
        z_pd32 = x_pd32(:,:,end);
    end 
    [breakpt_tmp32] = bfk_single_break_detect(N,T-est_break_pts(2),x_pd32,y_pd32,z_pd32);
    
    y_pd311_tmp = ymat(1:est_break_pts(1),:);
    x_pd311_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr311,bhat311] = bfk_post_break_est(N,est_break_pts(1),y_pd311_tmp,x_pd311_tmp);   
    y_pd312_tmp = ymat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp31,:);
    x_pd312_tmp = xmat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp31,:,:);
    [ssr312,bhat312] = bfk_post_break_est(N,breakpt_tmp31,y_pd312_tmp,x_pd312_tmp);
    y_pd313_tmp = ymat(est_break_pts(1)+breakpt_tmp31+1:est_break_pts(2),:);
    x_pd313_tmp = xmat(est_break_pts(1)+breakpt_tmp31+1:est_break_pts(2),:,:);
    [ssr313,bhat313] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1)-breakpt_tmp31,y_pd313_tmp,x_pd313_tmp);
    y_pd314_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd314_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr314,bhat314] = bfk_post_break_est(N,T-est_break_pts(2),y_pd314_tmp,x_pd314_tmp);
        
    y_pd321_tmp = ymat(1:est_break_pts(1),:);
    x_pd321_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr321,bhat321] = bfk_post_break_est(N,est_break_pts(1),y_pd321_tmp,x_pd321_tmp);   
    y_pd322_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd322_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr322,bhat322] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd322_tmp,x_pd322_tmp);
    y_pd323_tmp = ymat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:);
    x_pd323_tmp = xmat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:,:);
    [ssr323,bhat323] = bfk_post_break_est(N,breakpt_tmp32,y_pd323_tmp,x_pd323_tmp);
    y_pd324_tmp = ymat(est_break_pts(2)+breakpt_tmp32+1:T,:);
    x_pd324_tmp = xmat(est_break_pts(2)+breakpt_tmp32+1:T,:,:);
    [ssr324,bhat324] = bfk_post_break_est(N,T-est_break_pts(2)-breakpt_tmp32,y_pd324_tmp,x_pd324_tmp);

    ssr31 = ssr311 + ssr312 + ssr313 + ssr314;
    ssr32 = ssr321 + ssr322 + ssr323 + ssr324;    
    if ssr31 < ssr32
        est_break_pts(3) = est_break_pts(1)+breakpt_tmp31;
    else
        est_break_pts(3) = est_break_pts(2)+breakpt_tmp32;
    end   
    
elseif est_break_pts(2)-est_break_pts(1)<=2
    y_pd31 = ymat(1:est_break_pts(1),:);
    x_pd31 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd31 = x_pd31(:,:,1:end);
    else
        z_pd31 = x_pd31(:,:,end);
    end
    [breakpt_tmp31] = bfk_single_break_detect(N,est_break_pts(1),x_pd31,y_pd31,z_pd31);

    y_pd32 = ymat(est_break_pts(2)+1:T,:);
    x_pd32 = xmat(est_break_pts(2)+1:T,:,:);
    if ifull == 1
        z_pd32 = x_pd32(:,:,1:end);
    else
        z_pd32 = x_pd32(:,:,end);
    end
    [breakpt_tmp32] = bfk_single_break_detect(N,T-est_break_pts(2),x_pd32,y_pd32,z_pd32);
    
    y_pd311_tmp = ymat(1:breakpt_tmp31,:);
    x_pd311_tmp = xmat(1:breakpt_tmp31,:,:);
    [ssr311,bhat311] = bfk_post_break_est(N,breakpt_tmp31,y_pd311_tmp,x_pd311_tmp);   
    y_pd312_tmp = ymat(breakpt_tmp31+1:est_break_pts(1),:);
    x_pd312_tmp = xmat(breakpt_tmp31+1:est_break_pts(1),:,:);
    [ssr312,bhat312] = bfk_post_break_est(N,est_break_pts(1)-breakpt_tmp31,y_pd312_tmp,x_pd312_tmp);
    y_pd313_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd313_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr313,bhat313] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd313_tmp,x_pd313_tmp);
    y_pd314_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd314_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr314,bhat314] = bfk_post_break_est(N,T-est_break_pts(2),y_pd314_tmp,x_pd314_tmp);
        
    y_pd321_tmp = ymat(1:est_break_pts(1),:);
    x_pd321_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr321,bhat321] = bfk_post_break_est(N,est_break_pts(1),y_pd321_tmp,x_pd321_tmp);   
    y_pd322_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd322_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr322,bhat322] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd322_tmp,x_pd322_tmp);
    y_pd323_tmp = ymat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:);
    x_pd323_tmp = xmat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:,:);
    [ssr323,bhat323] = bfk_post_break_est(N,breakpt_tmp32,y_pd323_tmp,x_pd323_tmp);
    y_pd324_tmp = ymat(est_break_pts(2)+breakpt_tmp32+1:T,:);
    x_pd324_tmp = xmat(est_break_pts(2)+breakpt_tmp32+1:T,:,:);
    [ssr324,bhat324] = bfk_post_break_est(N,T-est_break_pts(2)-breakpt_tmp32,y_pd324_tmp,x_pd324_tmp);

    ssr31 = ssr311 + ssr312 + ssr313 + ssr314;
    ssr32 = ssr321 + ssr322 + ssr323 + ssr324;    
    if ssr31 < ssr32
        est_break_pts(3) = breakpt_tmp31;
    else
        est_break_pts(3) = est_break_pts(2)+breakpt_tmp32;
    end
    
    
elseif est_break_pts(2) >= T-2
    y_pd31 = ymat(1:est_break_pts(1),:);
    x_pd31 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd31 = x_pd31(:,:,1:end);
    else
        z_pd31 = x_pd31(:,:,end);
    end
    [breakpt_tmp31] = bfk_single_break_detect(N,est_break_pts(1),x_pd31,y_pd31,z_pd31);

    y_pd32 = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd32 = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    if ifull == 1
        z_pd32 = x_pd32(:,:,1:end);
    else
        z_pd32 = x_pd32(:,:,end);
    end
    [breakpt_tmp32] = bfk_single_break_detect(N,est_break_pts(2)-est_break_pts(1),x_pd32,y_pd32,z_pd32);
    
    y_pd311_tmp = ymat(1:breakpt_tmp31,:);
    x_pd311_tmp = xmat(1:breakpt_tmp31,:,:);
    [ssr311,bhat311] = bfk_post_break_est(N,breakpt_tmp31,y_pd311_tmp,x_pd311_tmp);   
    y_pd312_tmp = ymat(breakpt_tmp31+1:est_break_pts(1),:);
    x_pd312_tmp = xmat(breakpt_tmp31+1:est_break_pts(1),:,:);
    [ssr312,bhat312] = bfk_post_break_est(N,est_break_pts(1)-breakpt_tmp31,y_pd312_tmp,x_pd312_tmp);
    y_pd313_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd313_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr313,bhat313] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd313_tmp,x_pd313_tmp);
    y_pd314_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd314_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr314,bhat314] = bfk_post_break_est(N,T-est_break_pts(2),y_pd314_tmp,x_pd314_tmp);
        
    y_pd321_tmp = ymat(1:est_break_pts(1),:);
    x_pd321_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr321,bhat321] = bfk_post_break_est(N,est_break_pts(1),y_pd321_tmp,x_pd321_tmp);   
    y_pd322_tmp = ymat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp32,:);
    x_pd322_tmp = xmat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp32,:,:);
    [ssr322,bhat322] = bfk_post_break_est(N,breakpt_tmp32,y_pd322_tmp,x_pd322_tmp);
    y_pd323_tmp = ymat(est_break_pts(1)+breakpt_tmp32+1:est_break_pts(2),:);
    x_pd323_tmp = xmat(est_break_pts(1)+breakpt_tmp32+1:est_break_pts(2),:,:);
    [ssr323,bhat323] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1)-breakpt_tmp32,y_pd323_tmp,x_pd323_tmp);
    y_pd324_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd324_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr324,bhat324] = bfk_post_break_est(N,T-est_break_pts(2),y_pd324_tmp,x_pd324_tmp);

    ssr31 = ssr311 + ssr312 + ssr313 + ssr314;
    ssr32 = ssr321 + ssr322 + ssr323 + ssr324;    
    if ssr31 < ssr32
        est_break_pts(3) = breakpt_tmp31;
    else
        est_break_pts(3) = est_break_pts(1)+breakpt_tmp32;
    end   
    
else
    y_pd31 = ymat(1:est_break_pts(1),:);
    x_pd31 = xmat(1:est_break_pts(1),:,:);
    if ifull == 1
        z_pd31 = x_pd31(:,:,1:end);
    else
        z_pd31 = x_pd31(:,:,end);
    end
    [breakpt_tmp31] = bfk_single_break_detect(N,est_break_pts(1),x_pd31,y_pd31,z_pd31);

    y_pd32 = ymat(est_break_pts(2)+1:T,:);
    x_pd32 = xmat(est_break_pts(2)+1:T,:,:);
    if ifull == 1
        z_pd32 = x_pd32(:,:,1:end);
    else
        z_pd32 = x_pd32(:,:,end);
    end
    [breakpt_tmp32] = bfk_single_break_detect(N,T-est_break_pts(2),x_pd32,y_pd32,z_pd32);

    y_pd33 = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd33 = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    if ifull == 1
        z_pd33 = x_pd33(:,:,1:end);        
    else
        z_pd33 = x_pd33(:,:,end);        
    end    
    [breakpt_tmp33] = bfk_single_break_detect(N,est_break_pts(2)-est_break_pts(1),x_pd33,y_pd33,z_pd33);
      
    y_pd311_tmp = ymat(1:breakpt_tmp31,:);
    x_pd311_tmp = xmat(1:breakpt_tmp31,:,:);
    [ssr311,bhat311] = bfk_post_break_est(N,breakpt_tmp31,y_pd311_tmp,x_pd311_tmp);   
    y_pd312_tmp = ymat(breakpt_tmp31+1:est_break_pts(1),:);
    x_pd312_tmp = xmat(breakpt_tmp31+1:est_break_pts(1),:,:);
    [ssr312,bhat312] = bfk_post_break_est(N,est_break_pts(1)-breakpt_tmp31,y_pd312_tmp,x_pd312_tmp);
    y_pd313_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd313_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr313,bhat313] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd313_tmp,x_pd313_tmp);
    y_pd314_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd314_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr314,bhat314] = bfk_post_break_est(N,T-est_break_pts(2),y_pd314_tmp,x_pd314_tmp);
        
    y_pd321_tmp = ymat(1:est_break_pts(1),:);
    x_pd321_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr321,bhat321] = bfk_post_break_est(N,est_break_pts(1),y_pd321_tmp,x_pd321_tmp);   
    y_pd322_tmp = ymat(est_break_pts(1)+1:est_break_pts(2),:);
    x_pd322_tmp = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
    [ssr322,bhat322] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd322_tmp,x_pd322_tmp);
    y_pd323_tmp = ymat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:);
    x_pd323_tmp = xmat(est_break_pts(2)+1:est_break_pts(2)+breakpt_tmp32,:,:);
    [ssr323,bhat323] = bfk_post_break_est(N,breakpt_tmp32,y_pd323_tmp,x_pd323_tmp);
    y_pd324_tmp = ymat(est_break_pts(2)+breakpt_tmp32+1:T,:);
    x_pd324_tmp = xmat(est_break_pts(2)+breakpt_tmp32+1:T,:,:);
    [ssr324,bhat324] = bfk_post_break_est(N,T-est_break_pts(2)-breakpt_tmp32,y_pd324_tmp,x_pd324_tmp);
        
    y_pd331_tmp = ymat(1:est_break_pts(1),:);
    x_pd331_tmp = xmat(1:est_break_pts(1),:,:);
    [ssr331,bhat331] = bfk_post_break_est(N,est_break_pts(1),y_pd331_tmp,x_pd331_tmp);   
    y_pd332_tmp = ymat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp33,:);
    x_pd332_tmp = xmat(est_break_pts(1)+1:est_break_pts(1)+breakpt_tmp33,:,:);
    [ssr332,bhat332] = bfk_post_break_est(N,breakpt_tmp33,y_pd332_tmp,x_pd332_tmp);
    y_pd333_tmp = ymat(est_break_pts(1)+breakpt_tmp33+1:est_break_pts(2),:);
    x_pd333_tmp = xmat(est_break_pts(1)+breakpt_tmp33+1:est_break_pts(2),:,:);
    [ssr333,bhat333] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1)-breakpt_tmp33,y_pd333_tmp,x_pd333_tmp);
    y_pd334_tmp = ymat(est_break_pts(2)+1:T,:);
    x_pd334_tmp = xmat(est_break_pts(2)+1:T,:,:);
    [ssr334,bhat334] = bfk_post_break_est(N,T-est_break_pts(2),y_pd334_tmp,x_pd334_tmp);

    ssr31 = ssr311 + ssr312 + ssr313 + ssr314;
    ssr32 = ssr321 + ssr322 + ssr323 + ssr324;    
    ssr33 = ssr331 + ssr332 + ssr333 + ssr334; 
    [~, optbreak] = min([ssr31,ssr32,ssr33]);
    switch optbreak
        case 1
            est_break_pts(3) = breakpt_tmp31;
        case 2
            est_break_pts(3) = est_break_pts(2)+breakpt_tmp32;
        case 3
            est_break_pts(3) = est_break_pts(1)+breakpt_tmp33;
    end    
end



end
end


est_break_pts = sort(est_break_pts);

% y_pd1 = ymat(1:est_break_pts(1),:);
% x_pd1 = xmat(1:est_break_pts(1),:,:);
% [~,bhat1,sehat1] = bfk_post_break_est(N,est_break_pts(1),y_pd1,x_pd1);   
% y_pd2 = ymat(est_break_pts(1)+1:est_break_pts(2),:);
% x_pd2 = xmat(est_break_pts(1)+1:est_break_pts(2),:,:);
% [~,bhat2,sehat2] = bfk_post_break_est(N,est_break_pts(2)-est_break_pts(1),y_pd2,x_pd2);   
% y_pd3 = ymat(est_break_pts(2)+1:T,:);
% x_pd3 = xmat(est_break_pts(2)+1:T,:,:);
% [~,bhat3,sehat3] = bfk_post_break_est(N,T-est_break_pts(2),y_pd3,x_pd3);   
% 
% b1 = kron(bhat1,ones(1,est_break_pts(1)));
% b2 = kron(bhat2,ones(1,est_break_pts(2)-est_break_pts(1)));
% b3 = kron(bhat3,ones(1,T-est_break_pts(2)));
% se1 = kron(sehat1,ones(1,est_break_pts(1)));
% se2 = kron(sehat2,ones(1,est_break_pts(2)-est_break_pts(1)));
% se3 = kron(sehat3,ones(1,T-est_break_pts(2)));
% 
% bhat_tmp = reshape([b1 b2 b3]',T,K,N);
% bhat = reshape(permute(bhat_tmp,[1,3,2]),N*T,K);
% sehat_tmp = reshape([se1 se2 se3]',N*T,K);
% sehat = reshape(permute(sehat_tmp,[1,3,2]),N*T,K);

end