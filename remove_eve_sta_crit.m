function [ ] = remove_eve_sta_crit( infile1,outfile1,infile2,outfile2,crit,num_twin )
%remove_eve_sta_crit: remove the stations and sources whose records are 
%   smaller than crit. The removing process are iteratively run until 
%   all the unmatched stations and sources are removed. 
%
%   Input parameters:
%       infile:     Input file stores the original data
%       outfile:    Output file stores the corrected data
%       crit:       Critical Number 
%       num_twin:   Number of representive points of coda waves
%
%   Requied script: countmember(a,b)
%	Discription: Array a is a uniqued array of array b. Count number
%  		of each cell a in b
%
%   Written by Wei Wang @ SIO, 08/28/2017
%

crit=str2num(crit);
num_twin=str2num(num_twin);

% Read in the original file (infile) : Synthetic data part
fid_in = fopen(infile1, 'r');
C = textscan(fid_in, ['%d %f %f %f %f %f %f %s', repmat(' %f', 1, 2+num_twin)]);
fclose(fid_in);

numtol=length(C);  % total lines of input file
eveid_in=C{1,1};
for i=2:7
    data1(:,i-1)=(C{1,i});
end
stanm_in=C{1,8}; 
data2=C{1,9:numtol};
for i=9:numtol
    data2(:,i-8)=(C{1,i});
end

% Count the unique eveid number
eveid1=eveid_in; 
eveid_uniq1=unique(eveid1);
evenum1=countmember(eveid_uniq1,eveid1);

% Count the unique stanm number
stanm1=stanm_in;
stanm_uniq1=unique(stanm1);
stanum1=countmember(stanm_uniq1,stanm1);

% Remove all the stations and source number than crit
eve_nn1=find(evenum1<crit); 
sta_nn1=find(stanum1<crit);
indx=0;
while(length(eve_nn1>0) | length(sta_nn1>0))
% Deal with event first    
     % Find eveid (evenm<crit) in the eveid1
    [c,ia,ib]=intersect(eveid1,eveid_uniq1(eve_nn1));
    [tf,loc]=ismember(eveid1,eveid_uniq1(eve_nn1));
    idx=[1:length(eveid1)];
    idx=idx(tf);
    idx=idx(loc(tf));
    
    % Remove the corresponding line in eveid and stanm;
    eveid2=eveid1; eveid2(idx)=[];
    stanm2=stanm1; stanm2(idx)=[];
    
    eveid1=eveid2;
    stanm1=stanm2;
    
    % Count the unique stanm number
    eveid_uniq1=unique(eveid1);
    evenum1=countmember(eveid_uniq1,eveid1);

    % Count the unique stanm number
    stanm_uniq1=unique(stanm1);
    stanum1=countmember(stanm_uniq1,stanm1);
    
    eve_nn1=find(evenum1<crit); 
    sta_nn1=find(stanum1<crit);
    
% Deal with station second    
     % Find stanm (stanum<crit) in the stanm
    [c,ia,ib]=intersect(stanm1,stanm_uniq1(sta_nn1));
    [tf,loc]=ismember(stanm1,stanm_uniq1(sta_nn1));
    idx=[1:length(stanm1)];
    idx=idx(tf);
    idx=idx(loc(tf));
    
    % Remove the corresponding line in eveid and stanm;
    eveid2=eveid1; eveid2(idx)=[];
    stanm2=stanm1; stanm2(idx)=[];
    
    eveid1=eveid2;
    stanm1=stanm2;
    
    % Count the unique stanm number
    eveid_uniq1=unique(eveid1);
    evenum1=countmember(eveid_uniq1,eveid1);

    % Count the unique stanm number
    stanm_uniq1=unique(stanm1);
    stanum1=countmember(stanm_uniq1,stanm1);
    
    eve_nn1=find(evenum1<crit); 
    sta_nn1=find(stanum1<crit);    
    
    indx=indx+1
end

% Write out into new file
fid_out = fopen(outfile1, 'w');
for i=1:length(eveid1)
    nn=find(eveid_in==eveid1(i) & (strcmp(stanm_in,stanm1(i))==1));
    stanm_tmp=cell2mat(stanm_in(nn));
    fprintf(fid_out,'%10d %8.2f %7.2f %8.2f %6.2f %7.2f %4.2f %5s %7.2f %8.2f ',...
        eveid_in(nn),data1(nn,:),stanm_tmp,data2(nn,1:2));
    fprintf(fid_out,'%8.4f ',data2(nn,3:end));
    fprintf(fid_out,'\n');
end
fclose(fid_out);

%%
% Read in the original file (infile) : Observed data part
fid_in = fopen(infile2, 'r');
C = textscan(fid_in, ['%d %f %f %f %f %f %f %s', repmat(' %f', 1, 2+num_twin)]);
fclose(fid_in);

numtol=length(C);  % total lines of input file
eveid_in=C{1,1};
for i=2:7
    data1(:,i-1)=(C{1,i});
end
stanm_in=C{1,8}; 
data2=C{1,9:numtol};
for i=9:numtol
    data2(:,i-8)=(C{1,i});
end

% Write out into new file
fid_out = fopen(outfile2, 'w');
for i=1:length(eveid1)
    nn=find(eveid_in==eveid1(i) & (strcmp(stanm_in,stanm1(i))==1));
    stanm_tmp=cell2mat(stanm_in(nn));
    fprintf(fid_out,'%10d %8.2f %7.2f %8.2f %6.2f %7.2f %4.2f %5s %7.2f %8.2f ',...
        eveid_in(nn),data1(nn,:),stanm_tmp,data2(nn,1:2));
    fprintf(fid_out,'%8.4f ',data2(nn,3:end));
    fprintf(fid_out,'\n');
end
fclose(fid_out);


end

