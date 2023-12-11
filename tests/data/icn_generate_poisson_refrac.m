function T = icn_generate_poisson_refrac(interv,f,ref,opt)
% ICN_GENERATE_POISSON_REFRAC
% generate time series with a given refractory period
%
% T = icn_generate_poisson_refrac(interv,f,ref);
% T = icn_generate_poisson_refrac(interv,f,ref,'exact');
%
% inputs:
%    interv:   list of intervals (n-by-2 matrix)
%    f:        frequency in each interval (n-by-1 vector)
%    ref:      refractory period between consecutive spikes (seconds)
%
% outputs:
%    T:        time series such that the frequency in the interval
%              [interv(k,1), interv(k,2)] is f(k)
%              The option 'exact' ensures that there is exactly
%              (interv(k,2)-interv(k,1))*f(k) spikes in the interval

% Vincent Jacquemet
% Universite de Montreal

if nargin<4
    exact = 0;
else
    exact = ischar(opt);
end

T = [];
for k=1:size(interv,1)
    L = interv(k,2)-interv(k,1);
    if exact
        n = round(L*f(k));
        Ldim = L - (n-1)*ref;
        if Ldim>0
            fdim = n/Ldim;
            Tk = cumsum(-1/fdim*log(rand(n+1,1)));
            Tk = Tk/Tk(end)*Ldim + ref*(0:n)';
        else
            Tk = cumsum(-1/f(k)*log(rand(n+1,1)));
            Tk = Tk/Tk(end)*L;
        end
        Tk(end) = [];
    else
        n = max(1,round(1.5*L*f(k)));
        Tk = cumsum(ref-(1/f(k)-ref)*log(rand(n,1)));
        while Tk(end)<L
            Tk = [Tk; Tk(end)+cumsum(ref-(1/f(k)-ref)*log(rand(n,1)))];
        end
        Tk(Tk>=L) = [];
    end
    T = [T; interv(k,1) + Tk];
end

