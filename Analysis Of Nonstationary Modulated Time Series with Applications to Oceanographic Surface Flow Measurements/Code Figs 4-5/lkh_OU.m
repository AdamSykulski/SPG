function [out]= lkh_OU(x,xb,SX, ker, N,LB,UB,MF,ZEROF)
%This function computes the likelihood defined in the paper, using the
%kernel sent as the argument called ker. The latter needs to be precomputes
%via the call of the function coriolis_kernel(corfreqs). The reason why one
%cannot pass directly the coriolis frequencies to the likelihood function
%is that in the case where we want to minimize the said likelihood, the
%computation of the kernel is O(N^2) and does not depend on the parameters
%x. Since the rest of the likelihood calculation can be achieved in
%O(NlogN) we impose the separation of the kernel computation.
%
%Note: the coriolis frequencies are accounted for in the kernel k (which is
%complex-valued). But we still use 3 args when calling the function
%complexouacvs. The second arg is no longer the peak of the OU but a shift
%from the coriolis frequencies. For that reason the initial guess xb2(2) is
%to be set close to zero (no shift from the coriolis frequencies).
acv = complexouacvs(x(1:3) .* xb(1:3),N,1) .* ker;
ESF2=2*fft(acv)-acv(1); ESF3=abs(real(fftshift(ESF2)));
if ZEROF==0
    out=sum(log(ESF3([LB:MF-1 MF+1:UB])))+sum(SX([LB:MF-1 MF+1:UB])./ESF3([LB:MF-1 MF+1:UB]));
else
    out=sum(log(ESF3([LB:UB])))+sum(SX([LB:UB])./ESF3([LB:UB]));
end
end
