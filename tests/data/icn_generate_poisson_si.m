function T = icn_generate_poisson_si(interv,fT,SI,tau_s)
% ICN_GENERATE_POISSON_SI
% generate two time series with a given SI 
%
% T = icn_generate_poisson_si(interv,fT,SI,tau_s);
%
% inputs:
%    interv:    interval of the time series (seconds)
%    fT:        1-by-2 array containing the frequencies of the time series
%               T{1} and T{2} (Hz)
%    SI:        synchrony index the resulting T must have when T{2} is the
%               reference neuron
%    tau_s:     maximum delay for a coincidence (seconds)
%
% outputs:
%    T:         1-by-2 cell array containing the synchronized time series
%
% example:
%    T = icn_generate_poisson_si([0 200],[1 2],0.23,0.03);
%    SI2 = icn_synchrony_index(T,0.03);
%    SI2(1,2) == 0.23


if length(interv)<2
    interv = [0 interv];
end
if nargin<4
    tau_s = 0.04;
end

T = icn_timeseries_given_SI(SI,tau_s,fT,interv(2)-interv(1));
T{1} = T{1} + interv(1);
T{2} = T{2} + interv(1);



function T = icn_timeseries_given_SI(SI,tau_s,fT,duration)
% ICN_TIMESERIES_GIVEN_SI
% generate two time series with a given SI (T{2} is the reference)
%
% T = icn_timeseries_given_SI(SI,tau_s,fT,duration)
%
% inputs:
%    SI:        synchrony index the resulting T must have when T{2} is the
%               reference neuron
%    tau_s:     maximum delay for a coincidence (seconds)
%    fT:        1-by-2 array containing the frequencies of the time series
%               T{1} and T{2} (Hz)
%    duration:  length of the time series T{1} and T{2} (seconds)
%
% outputs:
%    T:         1-by-2 cell array containing the synchronized time series

T{1} = icn_generate_poisson_refrac([0 duration],fT(1),2*tau_s,'exact');

tau_j = 2*tau_s;
eps = 1e-9;

t = [T{1}-tau_s-tau_j;
            T{1}-tau_s-eps; T{1}-tau_s; T{1}-tau_s+eps; ...
            T{1}+tau_s-eps; T{1}+tau_s; T{1}+tau_s+eps; ...
            T{1}+tau_s+tau_j];
t = sort(t);

[S,p] = icn_synchrony_overlap_m(T{1}(:),t(:),tau_s,tau_j);
si=2*(S-p);

lim=1e-10;
beg0 = abs(si)<lim & [0;abs(si(1:end-1))]>lim   & [abs(si(2:end));inf]<lim; % pt=0 & bef~=0 & aft=0
end0 = abs(si)<lim & [inf;abs(si(1:end-1))]<lim & [abs(si(2:end));0]>lim;   % pt=0 & bef=0  & aft~=0
beg1 = si>(1-lim)  & [inf;si(1:end-1)]<(1-lim)  & [si(2:end);0]>(1-lim);    % pt=1 & bef<1  & aft=1
end1 = si>(1-lim)  & [0;si(1:end-1)]>(1-lim)    & [si(2:end);inf]<(1-lim);  % pt=1 & bef=1  & aft<1

interv0=[t(beg0) t(end0)];
interv0(any(interv0<0,2),:)=[];
interv0(any(interv0>duration,2),:)=[];
tot0=sum(interv0(:,2)-interv0(:,1));
tot0c=cumsum(interv0(:,2)-interv0(:,1));

interv1=[t(beg1) t(end1)];
interv1(any(interv1<0,2),:)=[];
interv1(any(interv1>duration,2),:)=[];
tot1=sum(interv1(:,2)-interv1(:,1));
tot1c=cumsum(interv1(:,2)-interv1(:,1));

n01=round(fT(2)*duration);
n1=round(SI*n01);
n0=n01-n1;

f0=n0/tot0;
f1=n1/tot1;

if isfinite(f0)
    Tt0 = icn_generate_poisson_refrac([0 tot0],f0,2*tau_s,'exact');
else
    Tt0=[];
end
if ~isempty(Tt0)
    [a,id] = sort([Tt0 ; tot0c]);
    q = (1:length(a))';
    q(id) = q;
    tc = cumsum(id>length(Tt0));        
    s = tc(q(1:length(Tt0)));
    iu = min(s+1,length(tot0c));
    adjust=[0;tot0c];
    Tt0=Tt0-adjust(iu)+interv0(iu,1);
end

if n1<=size(interv1,1)
    ints = randperm(size(interv1,1));
    ints = ints(1:n1);
    dt = interv1(ints,2)-interv1(ints,1);
    Tt1 = interv1(ints,1)+dt.*rand(n1,1);
else
%     disp('Not enough SI=1 intervals to place spikes')

    if isfinite(f1)
        Tt1 = icn_generate_poisson_refrac([0 tot1],f1,2*tau_s,'exact');
    %     Tt1 = icn_generate_poisson([0 tot1],f1,'exact');
    else
        Tt1=[];
    end
    if ~isempty(Tt1)
        [a,id] = sort([Tt1 ; tot1c]);
        q = (1:length(a))';
        q(id) = q;
        tc = cumsum(id>length(Tt1));        
        s = tc(q(1:length(Tt1)));
        iu = min(s+1,length(tot1c));
        adjust=[0;tot1c];
        Tt1=Tt1-adjust(iu)+interv1(iu,1);
    end   
end

T{2}=sort([Tt0;Tt1]);

