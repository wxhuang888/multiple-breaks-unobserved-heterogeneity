function EST_comp_bfk_Mknown(FILENAME, NEWFILENAME, M)


warning('off','all')
warning

%% load data
load(FILENAME);

%% sequential
CX = [X ones(N*T,1)]; 
K = size(CX,2);

xmat    = reshape(CX,T,N,K);
ymat    = reshape(y,T,N);
zmat    = xmat;
ifull   = 1;
% [est_break_pts, bhat, sehat] = bfk_multi_break_detect(N,T,ymat,xmat,zmat,ifull)

[est_break_pts] = Revise_bfk_multi_break_detect_M(N,T,ymat,xmat,zmat,ifull, M);
result.T0hat = est_break_pts;
result.T0 =regime0;
%% Save result 

save(NEWFILENAME,'result')    
    

end

