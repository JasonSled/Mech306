function DataMatrixOut = mech306_lab3a_1() 
    format long
    
    filename = 'Experiment-Numbers-Only.csv';
    fileID = fopen(filename);
    Data = textscan(fileID, '%s %s %s %s %s');
    fclose(fileID);
    
    timeStep = 2; %sec/measurement
    x = [3 54 104 155]; %mm
    x = .001*x; %m
    xtest1 = 1;
    Ts = 0; %C
    Ti = 25; %C
    
    %figure(1)
    %plot(Data(:, 1), Data(:, 2));
    
    DataMatrix = zeros(length(Data{1}), length(Data));
    for j = 1:length(Data)
        DataCol  = Data{j};
        for i = 1:length(Data{j})
            DataMatrix(i, j) = str2double(DataCol{i});
        end
    end
    
    DataMatrix(1, 1) = 0;
    DataMatrix(:, 1) = DataMatrix(:, 1).*timeStep;
    
    clf
    figure(1)
    hold
    %plot(DataMatrix(:, 1), DataMatrix(:, 2));
    %plot(DataMatrix(:, 1), DataMatrix(:, 3));
    %plot(DataMatrix(:, 1), DataMatrix(:, 4));
    %plot(DataMatrix(:, 1), DataMatrix(:, 5));
    plot(DataMatrix(:, 1), DataMatrix(:, xtest1+1));
    DataMatrixOut = DataMatrix;
    
    
    step = 0.00001;
    xtest = x(xtest1);
    
    %{
    TestAlphaAt500 = TestAlpha(xtest, Ts, Ti, step, 500);
    TestAlphaAt1000 = TestAlpha(xtest, Ts, Ti, step, 1000);
    TestAlphaAt1500 = TestAlpha(xtest, Ts, Ti, step, 1500);
    TestAlphaAt2000 = TestAlpha(xtest, Ts, Ti, step, 2000);
    
    
    scatter(TestAlphaAt500(:, 1), TestAlphaAt500(:, 2));
    scatter(TestAlphaAt1000(:, 1), TestAlphaAt1000(:, 2));
    scatter(TestAlphaAt1500(:, 1), TestAlphaAt1500(:, 2));
    scatter(TestAlphaAt2000(:, 1), TestAlphaAt2000(:, 2));
    
    %}
    
    alphaFinal = GetAlpha(0.00003,xtest, Ts, Ti, 0:2:2100);
    plot(alphaFinal(:, 1), alphaFinal(:, 2));
    
    title('Temperature over time at 3mm, alpha = 0.00003');
    ylabel('Temperature (C)');
    xlabel('Time (s)');
    
    %axis([0 2500 0 10])
    
end

function TestAlphaOut = TestAlpha(x, Ts, Ti, step, time)
    counts = 0.000001:step:0.0001;
    TestAlpha = ones(length(counts), 3);
    TestAlpha(:, 1) = TestAlpha(:, 1)*time;
    counter = 1;
    
    for k = counts
        T = (1 - erf(x/(2*sqrt(k*time))))*(Ts - Ti) + Ti;
        TestAlpha(counter, 2) = T;
        TestAlpha(counter, 3) = k;
        counter = counter + 1;
    end
    
    TestAlphaOut = TestAlpha;
end

function GetPlotForAlpha = GetAlpha(alpha, x, Ts, Ti, time)
    
    TestAlpha = ones(length(time), 1);
    counter = 1;
    
    for k = time
        T = (1 - erf(x/(2*sqrt(alpha*k))))*(Ts - Ti) + Ti;
        TestAlpha(counter, 1) = k;
        TestAlpha(counter, 2) = T;
        counter = counter + 1;
    end
    
    GetPlotForAlpha = TestAlpha;
end