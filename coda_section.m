%% Test new algorithm for coda_q 
% We first fit them applying linear regression and keep two points from the results.
clear all;
datadir='stp_data/';  % this sotres the sacfiles by the event id 
gridlist='evelist_test';  %  this is the event id list
stalist='sta_coda_q';	  % this is the prestored station name list
resdir='coda_step1_HZ_3Hz/'; % this is to sotre output coda infortatmion

% Filter parameter
[b,a]=butter(4,[0.04,0.08]);

% Window for taper (Unit: second)
len_taper=2.0;

% Parameter for twim window everage energy
ms_nt=512; nt_int=10;  % unit: point. this is for the coda smoothing (following Aki 1989)

% Parameter of time window for coda Q measurement
s_beg=50;  % (second) lapse time of coda window
s_twin=40; % (second) length of coda window


fid2=fopen(stalist,'r');
while(~feof(fid2))
    stanm=fgetl(fid2)
    tmpsta=strcat('*',stanm,'*HZ','.sac');
    stafile=strcat(resdir,stanm,'.txt');
    fid_sta=fopen(stafile,'w');
    fid1=fopen(gridlist,'r');
    while(~feof(fid1)) 
        tmpline=fgetl(fid1);
        griddir=strcat(datadir,tmpline);
        evelist=dir(griddir);
        evenum=length(evelist);
        for i=3:evenum
            evedir=strcat(griddir,'/',evelist(i).name,'/');
            sacnmlist=dir([evedir,tmpsta]);
            sacnum=length(sacnmlist);
            if sacnum ==1
                %arr_energy=zeros(s_nt,1);
                iflag=1;
                for isac=1:1
                    sacfile=strcat(evedir,sacnmlist(isac).name);
                    s=readsac(sacfile);
                    npts=s.NPTS;
                    if(npts<500 | isnan(npts))
                        iflag=0;
                        break;
                    end
                    delta=s.DELTA;
                    if(abs(delta-0.01)>1e-5 | abs(s.T1-5)>1e-2 | s.DIST>100)
                        iflag=0;
                        break
                    end
                    tt=s.B+(1:s.NPTS)*s.DELTA;
                    arr=s.DATA1; dtrend(arr);
                    [arr]=hann_taper(arr,delta,npts,len_taper);
                    arr_tmp = filter(b,a,arr);
                    nbeg1=floor(s.T1/delta+0.5)-200;nend1=nbeg1-1;
                    if(isnan(nbeg1))
                        iflag=0;
                        break;
                    end
                    noise_level=sum(arr_tmp(nbeg1:nend1).^2)/200;
                    tbeg=s.O+s_beg;  nbeg=floor(tbeg*100);
                    if(floor((tbeg+s_twin+10)*100)>npts)
                        iflag=0;
                        break
                    end
                    tend=nbeg+s_twin+5;
                    s_nt=floor(s_twin/nt_int/delta+0.5);
                    arr_energy=zeros(s_nt,1);
                    arr_ms=cal_ms_coda_q(arr_tmp,nbeg,s_nt,ms_nt,nt_int);
                    %[c_s,arrtmp]=check_after_shock(arr_ms,s_nt);

                    signal_level=sum(arr_ms)/s_nt;
                    snr=signal_level/noise_level;
                    if(snr<6) 
                        iflag=0;
                        break
                    end
                    arr_energy=arr_energy+arr_ms-noise_level;
                    [dist,az]=distance(s.EVLA,s.EVLO,s.STLA,s.STLO);
                    dist=s.DIST; az=s.AZ; mag=s.MAG;
                    if(isnan(dist))
                        iflag=0;
                    end
                end
                %tt2=50+(0:s_nt-1)*0.1;
                %plot(tt2,log10(arr_energy)); hold on;
                
                if(iflag==1 & arr_energy>0)
                    fprintf(fid_sta,'%s %8.2f %8.2f %8.2f %6.2f %7.2f %4.2f %s  %8.2f %8.2f ',...
                        evelist(i).name,s.EVLO,s.EVLA,s.EVDP,dist,az,mag,stanm,s.STLA,s.STLO);
                    fprintf(fid_sta,'%8.4f ',log10(arr_energy));
                    fprintf('%s  \n',evedir);
                    fprintf(fid_sta,'\n');
                end
            end
        end
    end
    fclose(fid_sta);
    fclose(fid1);
   
end
