function [post,classo] = post_est( N, a_hat, dat, this_group)
% INPUT:
%	T: number of time period
%	g_data: collect the data with group identity (group identity here is determined by CLasso)
%	class_a: coefficient estimated from CLasso

% OUTPUT: (each is a matrix containing the value of the coefficient, standard error and t-statistic.
%	classo: raw classo with bias correction
%	post: post lasso with bias correction


global p
index = 1:N;
g_index = index(this_group);
g_data = dat( ismember(dat.N, g_index), : );
%g_data_o = dato( ismember(dato.N, g_index), : );
classo.a= a_hat;
if all( logical( this_group ) == 0 )
    post_a = zeros(p, 1);
    post.a = post_a;
elseif sum(this_group) < p
    post_a = a_hat;
    post.a=post_a;
else
    [post_a,post_se,post_sig2]= b_est(g_data.y, g_data.X );
    post.a=post_a;
end
T = size(g_data,1) / size(g_index,2);
[vari_a] = var_PLS(T, a_hat', g_data.y, g_data.X  );
classo.std=sqrt(vari_a);


[vari_a_post] = var_PLS(T, post_a', g_data.y, g_data.X );

post.sig2 = post_sig2;
% post.std = post_se;
post.std = sqrt( diag( vari_a_post ) );

end
