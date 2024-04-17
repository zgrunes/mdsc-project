function importData_mlo(P, pth)
    % Read MLO monthly data into array.

    if nargin == 1
        pth = './data/raw';
    end

    dirIn = fullfile(pth, 'mlo');
    switch P.var
        case 'CO2'
            fileIn = 'co2_monthly.txt';
    end

    a = readcell(fullfile(dirIn, fileIn), 'NumHeaderLines', 159);

    nA = size(a, 1);
    dateStr = cell(nA, 1);
    for i = 1 : nA
        dateStr{i} = sprintf('%i%02i', a{i, 2}, a{i, 3});
    end
    idxT1 = find(strcmp(dateStr, P.tLim{1}), 1);
    idxT2 = find(strcmp(dateStr, P.tLim{2}), 1);
    x = cell2mat(a(idxT1 : idxT2, 11))';
    nD = 1;
    xystr = 'x1-1_y1-1'; % dummy string for lat/lon range

    dirOut = fullfile(pth, ...
                      'mlo', ...
                      'CO2', ...
                      [xystr '_' P.tLim{1} '-' P.tLim{2}] ...
                      );
    if ~isdir(dirOut)
        mkdir(dirOut)
    end
    save(fullfile(dirOut, 'dataGrid.mat'), 'nD')
    save(fullfile(dirOut, 'dataX.mat'), 'x')
end
