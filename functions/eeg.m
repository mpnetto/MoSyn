%% EEG
% 
%  Autor: Marcos Netto
%  Email: mpnetto88@gmail.com
%% Calcula a REDE EEG
% implementar comentario

%mex -v -largeArrayDims prototipo_mex.cpp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"


function Calc_Indexes
filename = 'C:\Users\Administrador\Documents\Marcos\Eeg\Teste\ANA PAULA MÃO 0097402.set.asc';
initialTime = 1;
endTime = 63000;
slidWindow =20;
tao = 3;
threshold = 0.8;
[data,labels] = readFile(filename);
numNodes = size(data,2);

global files;
guiHandle = app1;
guiFigureHandle  = guiHandle.Program;
uiwait(guiFigureHandle);
display(files);

filename = guiHandle.files(1);

tic;[data,labels] = readFile(filename);toc;
tic
[TVGArray,TVGArray_D,TVGWeig_Dir,motif] = prototipo_mex(data, initialTime, finalTime, numElec, slidWindow, tao, threshold);
toc
tic
ki_t = squeeze(mean(TVGArray,1));
k_t = sum(ki_t,2);

Degree_TVG= sum(TVGArray,2);
Degree_TVG = squeeze(Degree_TVG);
HubDeg = Degree_TVG > (mean (Degree_TVG, 1)+ 2*std(Degree_TVG,1));
                    
Degree_In = sum(TVGWeig_Dir,1);
Degree_In = squeeze(Degree_In);
HubDeg_In = Degree_In > (mean(Degree_In,1)+ 2*std(Degree_In,1));

Degree_Out = sum(TVGWeig_Dir,2);
Degree_Out = squeeze(Degree_Out);
HubDeg_Out = Degree_Out > (mean(Degree_Out,1)+ 2*std(Degree_Out,1));

REA = sum(TVGArray,3);
                
DegreeREA = sum(REA~=0,2);
HubSim = sum(HubDeg,2);
Hubs = HubSim./sum(HubSim,1);
HubIn = sum(HubDeg_In,2);
HubsIn = HubIn./sum(HubIn,1);
HubOut = sum(HubDeg_Out,2);
HubsOut = HubOut./sum(HubOut,1);

WeiDegREA = sum(sum(TVGArray,3),2);
WeiDegIn = sum(sum(TVGWeig_Dir,3),1)';
WeiDegOut = sum(sum(TVGWeig_Dir,3),2);

toc
tic
s = size(TVGArray_D);
s = s(3);

Cai_t = zeros(numElec,s,'single');
Li_t = zeros(numElec, numElec,s,'single');
T_t = zeros(numElec,s,'single');

for i = 1 : s
    a = TVGArray(:,:,i);
    Cai_t(:,i) = clustering_coef_bd(a);
    Li_t(:,:,i) = distance_bin(a);
    
    T_t(:,i) = sum(triu(a),2);
    
    %b = graph(a,'upper');
    %c = degree(b);
    %[d1, d2, d3]  = b.degree();
    %[f1, f2, f3] = b.eccentricity();
    %g = b.smallworldness();
    %h = b.pl();
    %i = b.diameter();
    %j = b.radius();
    %k = b.transitivity();
    %m = b.distance();
    %n = b.triangles();
end
Ca_t = mean(Cai_t,2); 
L_t = mean(Li_t,3);
toc

for i = 1:length(filename)

       MyListOfFiles{end+1,1} = filename(i);

end