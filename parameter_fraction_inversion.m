%% Do Coda Q inverstion and test the fraction of each parameters
% This part accounts for the importance of each part  
% 
clear;
filenm_syn='res_dir_HP_aplpha_1.5/freq_1.5/test_syn4_new';
filenm_dat='res_dir_HP_aplpha_1.5/freq_1.5/test_syn4_new';

sta_bad='sta_bad';
eve_bad='eve_bad';
syn_out='junk_syn';
dat_out1='junk_dat';
dat_out2='junk_dat2';

crit='10';
num_twin='400';

file_reorder_inversion(filenm_dat,dat_out2,eve_bad,sta_bad,'400','tmpeve','tmpsta');


%%
data=load(dat_out2);

% read in data 
eve=data(:,1);
eveid=data(:,3); evlo=data(:,4); evla=data(:,5); evdp=data(:,6);
sta=data(:,2); stla=data(:,10); stlo=data(:,11);
dist=data(:,7);
az=data(:,8);
mag=data(:,9);
amp=data(:,12:end);
clear data;

%%  This section is to use linear-regression to fit first half and second half to do coda-Q inv
tt1=50:0.1:69.9;    tt2=70:0.1:89.9;    tt3=50:0.1:89.9;
tt1_log=log10(tt1); tt2_log=log10(tt2); tt3_log=log10(tt3); 
dat1=amp(:,1:200); dat2=amp(:,201:400);
nrow=length(dat1); ncol=length(tt1);

alpha=1.5;

ncol=2;
b1_len=nrow*ncol; b1=zeros(b1_len,1);
b2_len=nrow*ncol; b2=zeros(b2_len,1);
b3_len=nrow*ncol; b3=zeros(b2_len,1);

tt1=50:0.1:69.9;    tt2=70:0.1:89.9;    tt3=50:0.1:89.9;
tt1_log=log10(tt1); tt2_log=log10(tt2); tt3_log=log10(tt3); 
dat1=amp(:,1:200); dat2=amp(:,201:400);
nrow=length(dat1); ncol=length(tt1);

ncol=2;
b1_len=nrow*ncol; b1=zeros(b1_len,1);
b2_len=nrow*ncol; b2=zeros(b2_len,1);
b3_len=nrow*ncol; b3=zeros(b2_len,1);
err1=zeros(nrow,1);  % error 
for i=1:nrow
    nbeg=(i-1)*ncol+1; nend=(i-1)*ncol+ncol; 
    
    tt=[50 70]; tt_log=log10(tt);
    A=zeros(400,2)+1; b=zeros(400,1);
    b=amp(i,:)'+alpha*tt3_log'; A(:,1)=tt3'*log10(exp(1));
    [x3,flag,err1(i,1)]=lsqr(A,b);
    b3(nbeg:nend)=x3(2)+x3(1)*tt'*log10(exp(1));
end
%%
tt=[50 70];
tt_log=log10(tt);

eve1=load('res_dir_HP_aplpha_1.5/freq_1.5/eve_res_both_3');
sta1=load('res_dir_HP_aplpha_1.5/freq_1.5/qc_res_both_3');

eve_res=eve1(:,6); qc_eve=eve1(:,7);
sta_res=sta1(:,5); qc_sta=sta1(:,4);

b_res=zeros(b3_len,1);  % best-fitting model
b_res1=zeros(b3_len,1);
b_res2=zeros(b3_len,1); b_res21=zeros(b3_len,1);
b_res3=zeros(b3_len,1); b_res31=zeros(b3_len,1);
b_res4=zeros(b3_len,1);
b_res5=zeros(b3_len,1);
b_res6=zeros(b3_len,1);
b_res7=zeros(b3_len,1);
b_res8=zeros(b3_len,1);
fid=fopen('test_bb1','w');
for i=1:nrow
    amp_tmp1=eve_res(eve(i))+sta_res(sta(i))-tt3*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)));
    amp_tmp2=amp(i,:)+alpha*tt3_log;
    err1(i,2)=norm(amp_tmp2-amp_tmp1)/norm(amp_tmp1);
    
    nbeg=(i-1)*ncol+1; nend=i*ncol;
    b_res(nbeg:nend)=b3(nbeg:nend);
    b_res1(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))-tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)));
    b_res2(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))...
                        -tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)))...
                        -alpha*tt_log';
    b_res21(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))...
                        -tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i))-mean(qc_eve)-mean(qc_sta));              
    b_res5(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))...
                        -tt'*log10(exp(1))*(qc_eve(eve(i))+mean(qc_sta));
    b_res4(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))...
                        -tt'*log10(exp(1))*(qc_sta(sta(i))+mean(qc_eve));
    b_res3(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i));
    b_res31(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))-tt'*log10(exp(1))*(mean(qc_sta)+mean(qc_eve));;
    b_res6(nbeg:nend)=mean(eve_res)+mean(sta_res)-tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)));
    b_res8(nbeg:nend)=eve_res(eve(i))+mean(sta_res)-tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)));
    b_res7(nbeg:nend)=mean(eve_res)+sta_res(sta(i))-tt'*log10(exp(1))*(qc_eve(eve(i))+qc_sta(sta(i)));
    fprintf(fid,'%5d %5d %5d %8.5f %8.5f %8.5f %8.5f\n',i,eve(i),sta(i),eve_res(eve(i)),qc_eve(eve(i)),sta_res(sta(i)),qc_sta(sta(i)));
end
fclose(fid)
aa=[b_res b_res1];
save('test_amp1','aa','-ascii');


err_tmp1=norm(b_res-b_res1)/norm(b_res);
err_tmp2=norm(b_res-b_res2)/norm(b_res); err_tmp21=norm(b_res-b_res21)/norm(b_res);
err_tmp3=norm(b_res-b_res3)/norm(b_res); err_tmp31=norm(b_res-b_res31)/norm(b_res);
err_tmp4=norm(b_res-b_res4)/norm(b_res);
err_tmp5=norm(b_res-b_res5)/norm(b_res);
err_tmp6=norm(b_res-b_res6)/norm(b_res);
err_tmp7=norm(b_res-b_res7)/norm(b_res);
err_tmp8=norm(b_res-b_res8)/norm(b_res);
%%
fprintf('%10.6f\n',err_tmp1);
fprintf('%10.6f\n',err_tmp2);
fprintf('%10.6f\n',err_tmp21);
fprintf('%10.6f\n',err_tmp3);
fprintf('%10.6f\n',err_tmp31);
fprintf('%10.6f\n',err_tmp4);
fprintf('%10.6f\n',err_tmp5);
fprintf('%10.6f\n',err_tmp6);
fprintf('%10.6f\n',err_tmp7);
fprintf('%10.6f\n',err_tmp8);
