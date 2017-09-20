function kt = createkt_asymptotes(ts,D,k,s,scalarlongtime,scalarshorttime)
% gives k(t), uses approximate k(t) and shot-time and long-time asymptotics
%scalarlong time ~10, scalarshortime ~0.001

kt = zeros(length(ts),1);
cutofflongtime = scalarlongtime*s*s/D;
cutoffshorttime = scalarshorttime*s*s/D;

for i = 1:length(ts)
    t = ts(i);
    
    if t<=cutoffshorttime
        
        % This is the short time asymptotic
        T = D.*t./s./s;
        
        if t == 0
            kt(i) = k;
        else
            if isinf(k)
                kt(i) = 2*pi*D*(1./(sqrt(pi*T)) + 1/2 - 1/4*sqrt(T/pi) + 1/8*T - 100/384*T.*sqrt(T/pi));
            else
                kt(i) = 2*pi*D*(1/2 - k/(2*pi*D)*sqrt(T/pi) - k/(4*pi*D)*T.*sqrt(T/pi) + k/(2*pi*D)*erfcx(k/(2*pi*D)*sqrt(T)).*(1 - pi*D/k + k/(2*pi*D)*T + 1/2*(k/(2*pi*D)*T).^2)  );
            end
        end
        
    elseif t>=cutofflongtime
        
        % This is the long time asymptotic
        
        B = 4*pi*D/k;
        g = -psi(1);
        
        C = D*t/s/s;
        % Eq longtassymptotic integral of Yogurtcu and Johnson 2015
        kt(i) = 4*pi*D*((1./(log(4*C)-2*g+B))-g*(1./(log(4*C)-2*g+B).^2)-1.331*(2./(log(4*C)-2*g+B).^3)+0.25*(2./(log(4*C)-2*g+B).^4));
        
    
    else
        kt(i)=-1;
    end
end
