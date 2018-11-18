function [ ] = file_reorder_inversion( infile,outfile,file_eve_bad,file_sta_bad,data_num,evefile,stafile )
%This scrips is used to reorder the original for inversion
%   All the staion names and eveids will be replaced by the index number  
%
%   Input parameters :
%       infile:     Input file stores the data from 'remove_eve_sta_crit'
%       outfile:    Output file stores reordered data for inversion
%       data_num:   Number of data point
%
%   Written by Wei Wang @ SIO, 09/08/2017
%


%infile='out_dat_test.txt';
%outfile='out_data_test_new.txt';
data_num=str2num(data_num);

% read in bad station
fid_sta = fopen(file_sta_bad,'r');
C1=textscan(fid_sta,'%s %f %f');
sta_bad=C1{1,1}; sta_bad_lon=C1{1,2}; sta_bad_lat=C1{1,3};
fclose(fid_sta);

% read in bad event
fid_eve = fopen(file_eve_bad,'r');
C1=textscan(fid_eve,'%d %f %f');
eve_bad=C1{1,1}; eve_bad_lon=C1{1,2}; eve_bad_lat=C1{1,3};
fclose(fid_eve);

% Read in the original file (infile)
fid_in = fopen(infile, 'r');
num=data_num+2;
C = textscan(fid_in, ['%d %f %f %f %f %f %f %s', repmat(' %f', 1, num)]);
fclose(fid_in);

numtol=length(C); 
eveid_in=C{1,1};
linetol=length(eveid_in);
for i=2:7
    data1(:,i-1)=(C{1,i});
end
stanm_in=C{1,8};
data2=C{1,9:12};
for i=9:numtol
    data2(:,i-8)=(C{1,i});
end

% remove bad staton and their corresponding lines
num_bad=length(sta_bad);
if(num_bad>0)
    for i=1:num_bad
        [sta_bad_log,sta_bad_indx]=ismember(stanm_in,sta_bad(i));
        %nn=find(sta_bad_indx==1);
        nn=find(sta_bad_indx==1 & data2(:,1)==sta_bad_lat(i) & data2(:,2)==sta_bad_lon(i));
        data1(nn,:)=[]; eveid_in(nn)=[];
        data2(nn,:)=[]; stanm_in(nn)=[];
    end
end
numtol=length(eveid_in);
linetol=length(eveid_in);

% remove bad event and their corresponding lines
num_bad=length(eve_bad);
if(num_bad>0)
    for i=1:num_bad
        [eve_bad_log,eve_bad_indx]=ismember(eveid_in,eve_bad(i));
        %nn=find(sta_bad_indx==1);
        nn=find(eve_bad_indx==1 & data1(:,2)==eve_bad_lat(i) & data1(:,1)==eve_bad_lon(i));
        data1(nn,:)=[]; eveid_in(nn)=[];
        data2(nn,:)=[]; stanm_in(nn)=[];
    end
end
numtol=length(eveid_in);
linetol=length(eveid_in);


% Reorder by the event id 
eveid1=(eveid_in);
stanm1=stanm_in(:);
data1=data1(:,:);
data2=data2(:,:);

eveid_uniq=unique(eveid1);
stanm_uniq=unique(stanm1);

[eveid_log,eveid_indx]=ismember(eveid1,eveid_uniq);
[stanm_log,stanm_indx]=ismember(stanm1,stanm_uniq);

% Write out the reordered file (outfile)
fid_out=fopen(outfile,'w');
for i=1:linetol
    fprintf(fid_out,'%4d %4d %8d %8.2f %7.2f %8.2f %6.2f %7.2f %4.2f %7.2f %8.2f ',...
        eveid_indx(i), stanm_indx(i), eveid1(i), data1(i,:),data2(i,1:2));
    fprintf(fid_out,'%8.4f ',data2(i,3:end));
    fprintf(fid_out,'\n');
end
fclose(fid_out);

% Write out the station name
fid_out=fopen(stafile,'w');
numsta=length(stanm_uniq)
for i=1:numsta
    tmpchar=cell2mat(stanm_uniq(i));
    fprintf(fid_out,'%5s\n',tmpchar);
end
fclose(fid_out);

% Write out the source id
fid_out=fopen(evefile,'w');
fprintf(fid_out,'%8d\n',eveid_uniq);
fclose(fid_out);
