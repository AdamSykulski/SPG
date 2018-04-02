function [out]=maternOUmodel(x,xb,SX,N,LB,UB,MF,ZEROF)
acv=maternacvs(x(4:6).*xb(4:6),N,1)+complexouacvs(x(1:3).*xb(1:3),N,1);
ESF2=2*fft(acv.*(1-(0:N-1)/N))-acv(1); ESF3=abs(real(fftshift(ESF2)));
if ZEROF==0
    out=sum(log(ESF3([LB:MF-1 MF+1:UB])))+sum(SX([LB:MF-1 MF+1:UB])./ESF3([LB:MF-1 MF+1:UB]));
else
    out=sum(log(ESF3([LB:UB])))+sum(SX([LB:UB])./ESF3([LB:UB]));
end