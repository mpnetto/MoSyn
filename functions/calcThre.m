%% Calculate Indexes
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementa    r comentario

%mex -v -largeArrayDims motif_sync_mex.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"

function threshold = calcThre(path, file,  wLength, wLag, taoMin,taoMax, mLag, porc)
    
    porc = porc/100;
    filename = char(strcat(path, file));

    [~,name,~] = fileparts(filename);

    [data,labels] = readFile(filename);
    
    initialTime = 0;
    
    endTime = size(data,1)-(wLength*2);

    numNodes = size(data,2);
    
    for i=1:numNodes
        newData(:,i)  = data(randperm(size(data,1)), i);
    end
    f = waitbar(0,'Please wait...');
    Q = calc_threshold_mex(newData, initialTime, endTime, numNodes, wLength, taoMin,taoMax);
    close(f);
    Sorted_Q = sort(Q(:));
    n = fix(length( Sorted_Q) *porc);
    
    el = Sorted_Q(n);
    
    threshold = double(el)  / double(wLength);
    disp(threshold);
    
end

