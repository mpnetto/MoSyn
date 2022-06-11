%% Calculate Indexes
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calculate 
% implementa    r comentario

%mex -v -largeArrayDims motif_sync_mex.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"

function [project] = calcIndexes(path, files, location, indexes, initialTime, finalTime, realTime, initialFrequency, finalFrequency, hasFrequency, removeTimeColumn, acqRate, slidWindow, taoMin,taoMax, threshold,tscLim,OutTvg)

    %multiWaitbar('Subjects',0, 'Color', 'g', 'CanCancel','on', 'CancelFcn', @(a,b) disp( ['Cancel ',a] ));
    
    numSubjects = length(files);
    
    subjects = zeros(1,numSubjects,'TVG');
    
    pFiles = strcat(path, files);
    
    for i = 1:numSubjects
        tic

        filename = pFiles{i};

        [~,name,ext] = fileparts(filename);

        if (ext == ".csv")
             [data,labels] = readCsvFile(filename);
        else
             [data,labels] = readFile(filename);
        end

        if(realTime)
            endTime = size(data,1)-slidWindow - taoMax -1;
        else
            endTime = finalTime;
        end
        
        if (endTime > (size(data,1)-slidWindow - taoMax -1))
            endTime = size(data,1)-slidWindow - taoMax -1;
            msg='Tamanho do arquivo menor do que o tempo final informado: Usando o tamanho como tempo final. Pressione OK para continuar. Arquivo:';
            msg=[msg filename];
            uiwait(msgbox(msg, 'Warning','modal'));
        end

        if(removeTimeColumn)
           data(:,1) = [];
           labels(1,:) = [];
        end

        if(hasFrequency)
            data = bandpass(data,[initialFrequency finalFrequency], acqRate);
        end

        numNodes = size(data,2);
        if(initialTime>0) 
            data=data(initialTime:end,:);
        end
        
       
          
        [TVGArray,TVGArray_D,TVGWeig_Dir, kkMotif] = motif_sync_mex(data, 0, endTime-initialTime, numNodes, slidWindow, taoMin,taoMax, threshold);
                 
        subjects(i) = TVG(name,path, TVGArray, TVGArray_D,TVGWeig_Dir, labels, indexes, initialTime, endTime, numNodes, slidWindow, taoMin,taoMax, threshold,tscLim,OutTvg,location);
            
        if subjects(i).runTVG()
            break
        end

        subjects(i).write();
        
         abort = multiWaitbar('Subjects', i/length(pFiles));
         
          if abort
              break
          end
        toc
    end

    project = Project(subjects, path, files, location, indexes);
    project.runProject();
    project.write();
   
    multiWaitbar('CloseAll');
end

