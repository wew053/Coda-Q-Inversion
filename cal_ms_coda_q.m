function [ arr_ms ] = cal_ms_coda_q( arrin,nbeg,s_nt,ms_nt,nt_int)
% cal_ms_coda_q Compute the mean-sqaure of time series
% arr_in: time series 
% nbeg: staring point of coda wave
% s_nt: time window length
% ms_nt: smoothing time window lenth
% nt_int: point skip for the smoothed time window

    arrtmp=arrin.^2;
    arr_ms=zeros(s_nt,1);
    for i=1:s_nt
        nbeg1=nbeg+(i-1)*nt_int-(ms_nt/2);
        nend1=nbeg1+(ms_nt)-1;
        arr_ms(i)=sum(arrtmp(nbeg1:nend1))/(2*floor(ms_nt/2));
    end
    arr_ms=(arr_ms);
end

