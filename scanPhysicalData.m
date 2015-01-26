%%  
function physicalDB = scanPhysicalData(filePath)

fid = fopen(filePath,'r');
tline = fgetl(fid);
lineCout = 0;
while tline(1) ~= -1
    lineCout = lineCout+1;
    dataStructure = textscan(tline,'%5s %1s %1s %f %f');
    tline = fgetl(fid);
end;
numberOfItems = lineCout;
frewind(fid);
physicalDB = struct;
physicalDB.filePath = filePath;
physicalDB.numberOfItems = numberOfItems;
tline = fgetl(fid);
height = zeros(numberOfItems,1);
weight = zeros(numberOfItems,1);
age = zeros(numberOfItems,1);
gender = zeros(numberOfItems,1);
lineCout = 0;
while tline(1) ~= -1
    lineCout = lineCout+1;
    dataStructure = textscan(tline,'%5s %1s %f %f %f');
    physicalDB.name{lineCout} = dataStructure{1};
    switch char(dataStructure{2})
        case 'F'
            gender(lineCout) = 0;
        otherwise
            gender(lineCout) = 1;
    end;
    age(lineCout) = dataStructure{3};
    height(lineCout) = dataStructure{4};
    weight(lineCout) = dataStructure{5};
    tline = fgetl(fid);
end;
physicalDB.gender = gender;
physicalDB.age = age;
physicalDB.height = height;
physicalDB.weight = weight;
fclose(fid);
