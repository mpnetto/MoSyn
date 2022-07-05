%% Calculate Indexes
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
  % This method takes a set of files, reads their data, transforms them into micro patterns (motifs), performs a synchrony analysis and transforms them into time-varying graphs.

%mex -v -largeArrayDims motif_sync_mex.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"

function [project] = calcIndexes(path, files, nodeLocation, indexes, initialTime, finalTime, realTime, initialFrequency, finalFrequency, hasFrequency, removeTimeColumn, acqRate, slidWindow, taoMin, taoMax, threshold, tscLim, calculateIndexes)

    project = -1;

    multiWaitbar('Subjects',0, 'Color', 'g', 'CanCancel','on', 'CancelFcn', @(a,b) disp( ['Cancel ',a] ));

    % The number of subjects is the number of files provided
    numSubjects = length(files);

    % Create a matrix with the (subject, value, index), in order to get the
    % value of each index selected for each subject
    subjects = zeros(1,numSubjects,'TVG');
    
    % Create an array with the name of the files
    pFiles = strcat(path, files);

    % Start timer to calculate the total calculation time
    totalTimeStart = tic;

    % Create an array with the name of the files
    for i = 1:numSubjects
       
        subjectTimeStart = tic;

        % Get the file name from the array
        filename = pFiles{i};

        % Extract the name and the extension
        [~,name,ext] = fileparts(filename);

        if (ext == ".mat")
             load (filename, 'TVGArray', 'TVGArray_D', 'TVGWeig_Dir','labels', 'initialTime', 'finalTime', 'realTime', 'endTime', 'initialFrequency', 'finalFrequency', 'hasFrequency', 'removeTimeColumn', 'acqRate', 'slidWindow', 'taoMin', 'taoMax', 'threshold', 'numNodes');

             if ~calculateIndexes
                 msgbox('asdasdasd', 'Warning','modal')
                 return
             end
        else

            % Check if the extension is from an csv file or not. Get the raw
            % data and labels from the file
            if (lower(ext) == ".csv")
                [data, labels] = readCsvFile(filename);
            elseif(lower(ext) == ".asc")
                 [data, labels] = readFile(filename);
            else
               return
            end
    
            % If the given start time is greater than 0, we remove from the data 
            % the rows corresponding to the time outside the range
            if(initialTime>0) 
                data=data(initialTime:end,:);
            end
    
            % the real end time is the numberof rows in the data, minus the slide window, minus the max TAO, minus 1
            realEndTime = size(data,1) - slidWindow - taoMax - 1;
    
            % Check if the the real time was provided, and that the real time is
            % less than size of the file. If not, the end time will be 
            if(realTime)
                endTime = realEndTime;
            elseif (realTime > realEndTime)
                endTime = realEndTime;
                msg=['Tamanho do arquivo menor do que o tempo final informado: Usando o tamanho como tempo final. Pressione OK para continuar. Arquivo:' filename];
                uiwait(msgbox(msg, 'Warning','modal'));
            else
                endTime = finalTime;
            end
    
           
             % Check if the "Allow time column" box was checked, and remove the first column
             % if so.
            if(removeTimeColumn)
               data(:,1) = [];
               labels(1,:) = [];
            end
    
            % Check if the "Allow Band Pass" box was checked and apply a
            % bandPass filter in the data
            if(hasFrequency)
                data = bandpass(data,[initialFrequency finalFrequency], acqRate);
            end
    
            % Number of nodes is equal to to the number of columns
            numNodes = size(data,2);
    
            % This method will transform the data into a motif matrix and calculate the TVG, the directed
            % TVG and the weighted directed TVG
            [TVGArray, TVGArray_D, TVGWeig_Dir, ~] = motif_sync_mex(data, 0, endTime-initialTime, numNodes, slidWindow, taoMin, taoMax, threshold);

            if (sum(TVGArray, 3) == 0)
                msgbox("Unable to sync data. Please try again with new parameters!")
                return
            end
    
            save ([path name '.tvg.mat'], 'TVGArray', 'TVGArray_D', 'TVGWeig_Dir', 'labels', 'initialTime', 'finalTime', 'realTime', 'endTime', 'initialFrequency', 'finalFrequency', 'hasFrequency', 'removeTimeColumn', 'acqRate', 'slidWindow', 'taoMin', 'taoMax', 'threshold', 'numNodes');
            
        end

        if calculateIndexes
            % This method will get the values data from the TVG and save into the
            % subject
            subjects(i) = TVG(name, path, TVGArray, TVGArray_D, TVGWeig_Dir, labels, indexes, initialTime, endTime, numNodes, slidWindow, taoMin, taoMax, threshold ,tscLim, nodeLocation);
                
            % This method will  caculate the indexes selected by the user in the UI
            subjects(i).runTVG()
            
            % After calculate the indexes theis method will write them into
            % files
            subjects(i).write();


            % Clar the TVgs for memory performance
            subjects(i).clearTVG();
               
        end

         % Finish timer and diplay calculate time of each subject
        subjectTimeEnd = toc(subjectTimeStart); 
        disp("Calculated time of subject " + name + " is: " + subjectTimeEnd + " seconds.");

        abort = multiWaitbar('Subjects', i/length(pFiles));

        if abort
            break
        end
        
    end

    if ~calculateIndexes
        msgbox({'The TVGs were generated and saved in .tvg.mat files in this path:'; path}, 'Warning','modal')
        totalTimeEnd = toc(totalTimeStart);
        disp("Total calculation time is: " + totalTimeEnd + " seconds.");
        return
    end

    %Configure and save project

    project = Project(subjects, path, files, nodeLocation, indexes);
    project.runProject();
    project.write();

     % Finish timer and diplays Total calculation time
    totalTimeEnd = toc(totalTimeStart);
    disp("Total calculation time is: " + totalTimeEnd + " seconds.");

    multiWaitbar('CloseAll');
   
end

