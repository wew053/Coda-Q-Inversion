%% Do Coda Q inverstion
% This part invert different part of coda-Q, sourre terms and station terms
% and results are stored int the dirctory

clear;
filenm_syn='res_dir_HP_aplpha_1.5/freq_3/test_syn4_new';
filenm_dat='res_dir_HP_aplpha_1.5/freq_3/test_dat4_new';

sta_bad='sta_bad';  	% the bad station 
eve_bad='eve_bad';  	% the eve band 
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
for i=1:nrow
    nbeg=(i-1)*ncol+1; nend=(i-1)*ncol+ncol; 
    
    tt=[50 70]; 
    tt_log=log10(tt);
    A=zeros(200,2)+1; b=zeros(200,1);
    b=dat1(i,:)'+alpha*tt1_log'; A(:,1)=tt1'*log10(exp(1));
    x1=lsqr(A,b);
    b1(nbeg:nend)=x1(2)+x1(1)*tt'*log10(exp(1))-alpha*tt_log';
    
    tt=[50 70];
    tt_log=log10(tt);
    A=zeros(200,2)+1; b=zeros(200,1);
    b=dat2(i,:)'+alpha*tt2_log'; A(:,1)=tt2'*log10(exp(1));
    x2=lsqr(A,b);
    b2(nbeg:nend)=x2(2)+x2(1)*tt'*log10(exp(1))-alpha*tt_log';
    
    
    tt=[50 70];
    A=zeros(400,2)+1; b=zeros(400,1);
    b=amp(i,:)'+alpha*tt3_log'; A(:,1)=tt3'*log10(exp(1));
    x3=lsqr(A,b);
    b3(nbeg:nend)=x3(2)+x3(1)*tt'*log10(exp(1))-alpha*tt_log';
end

%%  Invert the first part
tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+evenum+stanum;
b_nlen=b1_len;
val=1;
b=b1;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+eve(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    a4=sparse(nbeg:nend,evenum+stanum+evenum+sta(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    if(i==1)
        A=a1+a2+a3+a4;
    else
        A=A+a1+a2+a3+a4;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);


eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=x(evenum+stanum+1:evenum+stanum+evenum);
qc2_res=x(evenum+stanum+evenum+1:a_ncol);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_1'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_1'),'cc','-ascii');

%%  Invert the second part

tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+evenum+stanum;
b_nlen=b2_len;
val=1;
b=b2;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+eve(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    a4=sparse(nbeg:nend,evenum+stanum+evenum+sta(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    if(i==1)
        A=a1+a2+a3+a4;
    else
        A=A+a1+a2+a3+a4;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);


eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=x(evenum+stanum+1:evenum+stanum+evenum);
qc2_res=x(evenum+stanum+evenum+1:a_ncol);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_2'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_2'),'cc','-ascii');


%%  Invert the whole part

tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+evenum+stanum;
b_nlen=b3_len;
val=1;
b=b3;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+eve(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    a4=sparse(nbeg:nend,evenum+stanum+evenum+sta(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    if(i==1)
        A=a1+a2+a3+a4;
    else
        A=A+a1+a2+a3+a4;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);


eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=x(evenum+stanum+1:evenum+stanum+evenum);
qc2_res=x(evenum+stanum+evenum+1:a_ncol);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_3'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_3'),'cc','-ascii');
%%
b_res1=zeros(b_nlen,1);  % best-fitting model
b_res2=zeros(b_nlen,1);  % best-fitting model
tt=[50 70];
tt_log=log10(tt);
fid=fopen('test_bb','w');
for i=1:nrow

    nbeg=1+(i-1)*ncol; nend=i*ncol;
    b_res1(nbeg:nend)=b3(nbeg:nend)+alpha*tt_log';
    b_res2(nbeg:nend)=eve_res(eve(i))+sta_res(sta(i))-tt'*log10(exp(1))*(qc1_res(eve(i))+qc2_res(sta(i)));
    
    fprintf(fid,'%5d %5d %5d %8.5f %8.5f %8.5f %8.5f\n',i,eve(i),sta(i),eve_res(eve(i)),qc1_res(eve(i)),sta_res(sta(i)),qc2_res(sta(i)));
end
fclose(fid);
err_tmp1=norm(b_res1-b_res2)/norm(b_res1);
aa=[b_res1 b_res2];
save('test_amp','aa','-ascii');


%% invert for all station side coda-Q
tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+stanum;
b_nlen=b3_len;
val=1;
b=b3;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+sta(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    if(i==1)
        A=a1+a2+a3;
    else
        A=A+a1+a2+a3;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
save('test_sta','aa','-ascii');
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);

eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=zeros(evenum,1);
qc2_res=x(evenum+stanum+1:a_ncol);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_41'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_41'),'cc','-ascii');


%% invert for  source side coda-Q
tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+evenum;
b_nlen=b3_len;
val=1;
b=b3;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+eve(i),-tt*log10(exp(1)),b_nlen,a_ncol);
    if(i==1)
        A=a1+a2+a3;
    else
        A=A+a1+a2+a3;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);


eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=x(evenum+stanum+1:evenum+stanum+evenum);
qc2_res=zeros(stanum,1);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_51'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_51'),'cc','-ascii');


%% invert for one coda-Q
tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum+1;
b_nlen=b3_len;
val=1;
b=b3;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    a3=sparse(nbeg:nend,evenum+stanum+1,-tt*log10(exp(1)),b_nlen,a_ncol);

    if(i==1)
        A=a1+a2+a3;
    else
        A=A+a1+a2+a3;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);

eve_res=x(1:evenum);
sta_res=x(evenum+1:evenum+stanum);
qc1_res=x(evenum+stanum+1:evenum+stanum+1);
qc2_res=zeros(stanum,1);

% Write out results
% write out station terms and coda Q
qc=1./(qc2_res/2/pi/3);
aa=zeros(nrow,5);
aa(:,1)=sta; aa(:,2)=stlo; aa(:,3)=stla; 
aa(:,4)=qc2_res(sta); aa(:,5)=sta_res(sta);
bb=unique(aa,'rows');
cc=zeros(length(bb),5);
cc(:,1:5)=bb; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/qc_res_both_6'),'cc','-ascii');
% Write out error results
aa=zeros(nrow,5);
aa(:,1)=eveid; aa(:,2)=evlo; aa(:,3)=evla; aa(:,4)=evdp;
aa(:,5)=mag; 
bb=unique(aa,'rows');
cc=zeros(evenum,7);
cc(:,1:5)=bb; cc(:,6)=eve_res;cc(:,7)=qc1_res; 
save(fullfile(pwd,'res_dir_HP_aplpha_1.5/freq_3/eve_res_both_6'),'cc','-ascii');



%% invert for none coda-Q
tt=[50 70];
tt_log=log10(tt);

evenum=max(eve);
stanum=max(sta);
a_ncol=evenum+stanum;
b_nlen=b3_len;
val=1;
b=b3;
A=zeros(b_nlen,a_ncol);
for i=1:nrow
    nbeg=1+(i-1)*ncol;
    nend=i*ncol;
    b(nbeg:nend)=b(nbeg:nend)+alpha*tt_log';
    a1=sparse(nbeg:nend,eve(i),val,b_nlen,a_ncol);
    a2=sparse(nbeg:nend,evenum+sta(i),val,b_nlen,a_ncol);
    
    if(i==1)
        A=a1+a2;
    else
        A=A+a1+a2;
    end
end

% invert result
[x,flag,res]=lsqr(A,b,1e-8,2000);
res1=res*norm(b);
fprintf('%10.6f %12.6f\n',res,res1);
