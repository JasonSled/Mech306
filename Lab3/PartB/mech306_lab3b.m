function answer = mech306_lab3b()
    clear
    clf
    
    l1= .158; %black rod, len, m
    d1 = .009; %diam, m
    l2 = .157; %silver rod, len, m
    d2 = .009; %m
    
    %black paint, low fan
    v = 2.16; %m/s
    T = [39.458 42.220 49.531 60.179]; %C
    T = T+273.15;
    Tsurr = 25+273.15; %C
    V = [.99 4 7 10.01]; %V
    i = [.05 .26 .46 .65]; %A
    t = [60 60]; %sec
    x = .158; %distance from fan, m
    A = d1*l1*pi+pi*(d1^2)/2; %surface area of rod, m^2
    
    %{
    Since dT/dt = 0, we can simplify our equation to:
    q{elec} = q{radiation} + q{convection}
    (I)(V) = (epsilon)(sigma)(A)(T^4 - Tsurr^4) + (h)(A)(T - Tsurr)
    this can be simplified to:
    (I)(V) = (a)(T^4) - (a)(Tsurr^4) + (b)(T) - (b)(Tsurr)
    %}
    
    %func getCoeff solves for a & b, given a set of IV, a set of T, and Tsurr
    %and plot theoretical values in blue, real value in red
    [a, b] = getCoeff(T, V, i, Tsurr);
    
    %NOTE that our y = IV, x = T is very linear, 
    %so trying to make it fit a n=4 polynomial is either going to give you:
    %1. weird a & b
    %2. weird theoretical plot (r value is low)
    
    %NOTE that the way matlab interpolate plots, there can be multiple values
    %run this code multiple times and compare h to h from Whitaker approximation
    
    %func solveForEpiAndHFromCoeff solve (a) = (epsilon)(sigma)(A) and (b) = (h)(A)
    %for epsilon and h, respectively
    answer = solveForEpiAndHFromCoeff(a, b, A);
end

function [coeff1, coeff2] = getCoeff(T, V, i, Tsurr)
    
    x = T';
    y = i.*V;
    y = y';
    
    ft = fittype( @(a,b,c,d,x) a*x.^4+b*x+a*c+b*d, 'problem', {'c', 'd'});
    [fitobj,~,~,~]=fit(x,y,ft, 'problem', {-(Tsurr^4), -Tsurr}); %[fitobj, goodness, output, convmsg]=fit(x,y,ft)
     coeff1=fitobj.a;
      coeff2=fitobj.b;
      coeff3=fitobj.c;
      coeff4=fitobj.d;

    %{
    p = polyfit(x, y, 4);
    coeff1 = p(1);
    coeff2 = p(3);
    coeff3 = p(5);
    %}
     yfit=coeff1*x.^4+coeff2*x-coeff1*(Tsurr^4)-coeff2*Tsurr;     %Remember dot notation
     
     clf
     figure(1)
     hold
     plot(x, y, '-r');
     plot(x, yfit, '-b');
end

function answer = solveForEpiAndHFromCoeff(a, b, A)
    answer = zeros(2);
    answer(1) = a/(5.67*10^(-8)*A);
    answer(2) = b/A;
end
