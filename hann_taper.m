function [ arrout ] = hann_taper( arrin,delta,npts,len)
%hann_taper taper the time series with time windon 'len second'
    nlen=floor(len/delta);
    arg=(0:nlen-1)*pi/2/nlen;
    tt=sin(arg);
    arrout=arrin;
    arrout(1:nlen)=arrout(1:nlen).*tt';
    tt=flip(tt);
    arrout(npts-nlen+1:npts)=arrout(npts-nlen+1:npts).*tt';

end

