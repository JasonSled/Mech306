function R = mech306_lab3a_2()
    format long
    
    filename = 'Experiment-Numbers-Only.csv';
    fileID = fopen(filename);
    Data = textscan(fileID, '%s %s %s %s %s');
    fclose(fileID);
    
    timeStep = 2; %sec/measurement
    x = [3 54 104 155]; %mm
    x = .001*x; %m
    Tsurr = 273.15+25; %C
    m = .215; %kg
    C = 4.12 * (10^(-4)); %J/(kg*K)
    
    DataMatrix = zeros(length(Data{1}), length(Data));
    for j = 1:length(Data)
        DataCol  = Data{j};
        for i = 1:length(Data{j})
            DataMatrix(i, j) = str2double(DataCol{i});
        end
    end
    
    DataMatrix(1, 1) = 0;
    %DataMatrixtemp = DataMatrix(1:700, :);
    DataMatrixtemp = [DataMatrix(1:100, :); DataMatrix(150:700, :)];
    DataMatrix = DataMatrixtemp;
    DataMatrix(:, 1) = DataMatrix(:, 1).*timeStep;
    
    avgDataMatrix = zeros(size(DataMatrix));
    avgDataMatrix(:, 1) = DataMatrix(:, 1);
    avgDataMatrix(:, 2) = (DataMatrix(:, 2)+DataMatrix(:, 3)+DataMatrix(:, 4)+DataMatrix(:, 5))./4;
    avgDataMatrix(:, 2) = avgDataMatrix(:, 2)+273.15; %C to K
    
    step = 10;
    
    Rtotal = 0;
    counterAll = 0;
    counterGood = 0;
    RKeepTrackAll = [];
    RKeepTrackGood = [];
    for k = 2:step:length(avgDataMatrix(:, 2))-1
        frontdT = ((avgDataMatrix(k+1, 2) - avgDataMatrix(k, 2))/(avgDataMatrix(k+1, 1) - avgDataMatrix(k, 1)));
        backdT = ((avgDataMatrix(k, 2) - avgDataMatrix(k-1, 2))/(avgDataMatrix(k, 1) - avgDataMatrix(k-1, 1)));
        avgdT = (frontdT + backdT) / 2;
        Rtemp = (Tsurr - avgDataMatrix(k, 2)) / (m*C*avgdT);
        RKeepTrackAll(counterAll+1) = Rtemp;
        counterAll = counterAll + 1;
        if (Rtemp < 0 && Rtemp ~= inf && Rtemp ~= -inf)
            disp(Rtemp);
            RKeepTrackGood(counterGood+1) = Rtemp;
            Rtotal = Rtotal + Rtemp;
            counterGood = counterGood + 1;
        end
    end
    %R = Rtotal/counterGood;
    
    RremoveRows = RKeepTrackGood;
    RremoveRows(38) = [];
    RremoveRows(29) = [];
    RremoveRows(14) = [];
    
    R = sum(RremoveRows)/(counterGood - 3);
    
    clf
    figure(1)
    hold
    plot(avgDataMatrix(:, 1), avgDataMatrix(:, 2));
    
    %clf
    figure(2)
    hold
    plot(RKeepTrackAll);
    plot(RKeepTrackGood);
    plot(RremoveRows);
    yline(R);
end