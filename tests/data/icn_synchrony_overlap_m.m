function [S,p,a,b] = icn_synchrony_overlap_m(x,y,tau_s,tau_j)

% typically, length(x) > length(y)
% tau_s -> x; tau_j -> y

n = length(y);

% determine the list of disjoint intervals (x)
a = x(1)-tau_s;
b = x(1)+tau_s;
for i=2:length(x)
    if b(end)<x(i)-tau_s
        a = [a; x(i)-tau_s];
        b = [b; x(i)+tau_s];
    else
        b(end) = x(i)+tau_s;
    end
end

% compute intersections
p = zeros(n,1);
S = zeros(n,1);
j0 = 1;
for i=1:n
    c = [];
    d = [];
    while j0<length(a) && a(j0+1)<y(i)
        j0 = j0 + 1;
    end
    if a(j0)<= y(i) && y(i)<=b(j0)
        S(i) = 1;
    end
    for j=j0:length(a)
        if a(j)>=y(i)+tau_j
            break
        end
        cj = max(y(i)-tau_j,a(j));
        dj = min(y(i)+tau_j,b(j));
        if cj<dj
            c = [c; cj];
            d = [d; dj];
        end
    end
    for j=j0-1:-1:1
        if b(j)<=y(i)-tau_j
            break
        end
        cj = max(y(i)-tau_j,a(j));
        dj = min(y(i)+tau_j,b(j));
        if cj<dj
            c = [c; cj];
            d = [d; dj];
        end
    end
    if ~isempty(c)
        p(i) = sum(d-c)/2/tau_j;
    end
end

