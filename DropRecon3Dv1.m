classdef DropRecon3Dv1
    %DROPRECON3DV1 contains functions for constructing and analysing a
    % point cloud from a 3D tif-stack of a surface labeled drop or similar
    % object.
    %
    % Type 'DropRecon3Dv1.run()' to run.
    %
    % Copyright (c) 2016,	Elijah Shelton & Friedhelm Serwane
    %                       University of California, Santa Barbara
    % Last Modified July 2017
    %
    % License:              For Academic Use Only
    % Contact:              https://github.com/ElijahShelton
    % Coded by:             Elijah Shelton & Friedhelm Serwane
    
    properties
    end
    
    methods (Static)
        function [] = readme()
            readme_str = ['DropRecon Version 1.0 // 17 July 2017 \n',...
                '\n',...
                'GENERAL USAGE NOTES\n',...
                '---------------------------------------------------------------\n',...
                'This program has been developed for (1) constructing \n',...
                'point cloud models from tif-stacks of surface labeled \n',...
                'objects, such as deformable micro-droplets, and (2) \n',...
                'measuring the mean and gaussian curvature of their surfaces.\n',...
                '\n',...
                '\n',...
                'QUICK START:\n',...
                '\n',...
                '1) In order to run drop-recon, you must first compile mex \n',...
                'files from the following cpp files:\n',...
                '\n',...
                '\tfitGauss.cpp \t\t\t (required) \n',...
                '\tfitPolu1P10P01P00.cpp \t\t (required) \n',...
                '\tfitPoly2.cpp \t\t\t (required) \n',...
                '\tfitPoly2X.cpp \t\t\t (required) \n',...
                '\t- - - - - - - - - - - - - - - - - - - - - - \n',...
                '\tsteerableDetectior3D.cpp \t (optional) \n',...
                '\n',...
                '2) In MATLAB, cd into directory containing DropRecon3Dv1.m.\n',...
                '\n',...
                '3) To run, type ''DropRecon3Dv1.run'' and press enter.\n',...
                '\n',...
                'NOTE: MEX must be configured to compile C++ files.\n',...
                'In the matlab command line, type: mex -setup cpp\n\n',...
                'NOTE: If the required mex files are not present in the \n',...
                'local directory, the program will attempt to compile them \n',...
                'from the cpp files in ./lib/. After the mex files have \n',...
                'been compiled, you no longer need ./lib/ to run the program.\n',...
                '\n',...
                '4) You will be asked to locate the directory containing the \n',...
                'tiff stack or tiff stacks to be analyzed.\n',...
                '\n',...
                'DEMO: Select the directory ./tifs for test/ as a demo).\n',...
                '\n',...
                '5) You will be prompted to provide the voxel sizes in microns.\n',...
                '\n',...
                'DEMO: See /drop_info.txt for the voxel sizes for the demo drops\n',...
                '\n',...
                'NOTE: Tifs must have the same set of voxel sizes in order to \n',...
                'analyze them together in the same directory.\n',...
                '\n',...
                '6) You will be asked to select either CUSTOM or DEFAULT \n',...
                'parameters. DEFAULT is recommended.\n',...
                '\n',...
                'To find this information in matlab, type ''DropRecon3Dv1.readme'' \n',...
                'and press enter.\n',...
                '---------------------------------------------------------------\n',...
                'Please direct questions, comments, and suggestions to Elijah at\n',...
                'https://github.com/ElijahShelton',...
                '\n',...
                'This code was developed by Elijah Shelton and Friedhelm Serwane \n',...
                'at UC Santa Barbara, with input from Otger Campas.\n',...
                'Elijah Shelton received support from NSF GRFP while developing \n',...
                'this work.\n'];
            clc
            fprintf(readme_str)
            
        end
        function [] = turnOffWarnings()
            
            id = {'MATLAB:mir_warning_maybe_uninitialized_temporary';...
                'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning';...
                'MATLAB:imagesci:tiffmexutils:libtiffWarning'};
            numWarn = length(id);
            for n = 1:numWarn
                warning('off',id{n});
            end
            
        end
        function [] = compile_fitGauss_pc()
            if ispc % If running on a pc
                lsp = 'C:\gsl\lib';
                hsp = 'C:\gsl\include';
                if ~exist(lsp,'dir') || ~exist(hsp,'dir')
                    warning('You appear to not have the GSL libraries required for compilation.')
                end
                LDFLAGS = ['-L',lsp];
                CPPFLAGS = ['-I',hsp];
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitGauss.cpp');
                fprintf('\n''fitGauss.cpp'' Compilation Successful!\n')
            else
                fprintf(['This function is for compiling with a pc.\n',...
                    'MEX files were NOT created.'])
            end
        end
        function [] = compile_fitPoly_pc()
            if ispc % If running on a pc
                lsp = 'C:\gsl\lib';
                hsp = 'C:\gsl\include';
                if ~exist(lsp,'dir') || ~exist(hsp,'dir')
                    warning('You appear to not have the GSL libraries required for compilation.')
                end
                LDFLAGS = '-L/usr/local/opt/gsl@1/lib/';
                CPPFLAGS = '-I/usr/local/opt/gsl@1/include/';
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly1P10P01P00.cpp');
                fprintf('\n''fitPoly1P10P01P00.cpp'' Compilation Successful!\n')
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly2X.cpp');
                fprintf('\n''fitPoly2X.cpp'' Compilation Successful!\n')
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly2.cpp');
                fprintf('\n''fitPoly2.cpp'' Compilation Successful!\n')
                
            else
                fprintf(['This function is for compiling with a mac.\n',...
                    'MEX files were NOT created.'])
            end
        end
        function [] = compile_fitGauss_mac()
            if ismac % If running on a mac
                lsp = '/usr/local/opt/gsl@1/lib/';
                hsp = '/usr/local/opt/gsl@1/include/';
                if ~exist(lsp,'dir') || ~exist(hsp,'dir')
                    warning('You appear to not have the GSL libraries required for compilation.')
                    fprintf('Take the following steps to fix this problem:\n')
                    fprintf('1) Install Homebrew ( available at http://brew.sh ).\n')
                    fprintf('2) In the terminal, type: brew search gsl \n')
                    fprintf('3) In the terminal, type: brew install gsl@1 \n')
                    fprintf('4) In the matlab prompt, type: DropRecon3Dv1.compile_fitGauss_mac \n')
                end
                LDFLAGS = '-L/usr/local/opt/gsl@1/lib/';
                CPPFLAGS = '-I/usr/local/opt/gsl@1/include/';
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitGauss.cpp');
                fprintf('\n''fitGauss.cpp'' Compilation Successful!\n')
            else
                fprintf(['This function is for compiling with a mac.\n',...
                    'MEX files were NOT created.'])
            end
        end
        function [] = compile_fitPoly_mac()
            if ismac % If running on a mac
                lsp = '/usr/local/opt/gsl@1/lib/';
                hsp = '/usr/local/opt/gsl@1/include/';
                if ~exist(lsp,'dir') || ~exist(hsp,'dir')
                    warning('You appear to not have the GSL libraries required for compilation.')
                    fprintf('Take the following steps to fix this problem:\n')
                    fprintf('1) Install Homebrew ( available at http://brew.sh ).\n')
                    fprintf('2) In the terminal, type: brew search gsl \n')
                    fprintf('3) In the terminal, type: brew install gsl@1 \n')
                    fprintf('4) In the matlab prompt, type: DropRecon3Dv1.compile_fitPoly_mac \n')
                end
                LDFLAGS = '-L/usr/local/opt/gsl@1/lib/';
                CPPFLAGS = '-I/usr/local/opt/gsl@1/include/';
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly1P10P01P00.cpp');
                fprintf('\n''fitPoly1P10P01P00.cpp'' Compilation Successful!\n')
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly2X.cpp');
                fprintf('\n''fitPoly2X.cpp'' Compilation Successful!\n')
                mex(LDFLAGS,CPPFLAGS,'-lgslcblas','-lgsl','lib/fitPoly2.cpp');
                fprintf('\n''fitPoly2.cpp'' Compilation Successful!\n')
                
            else
                fprintf(['This function is for compiling with a mac.\n',...
                    'MEX files were NOT created.'])
            end
        end
        function [] = compileMex()
            if ismac
                DropRecon3Dv1.compile_fitGauss_mac();
                DropRecon3Dv1.compile_fitPoly_mac();
            elseif ispc
                DropRecon3Dv1.compile_fitGauss_pc();
                DropRecon3Dv1.compile_fitPoly_pc();
            end
        end
        function [] = mexCheck()
            DropRecon3Dv1.turnOffWarnings();
            
            if ~exist('fitPoly2','file')||~exist('fitPoly2X','file')||...
                    ~exist('fitPoly1P10P01P00','file')||...
                    ~exist('fitGauss','file')
                DropRecon3Dv1.compileMex();
            end
        end
        function [drop] = run()
            % DROP = RUN()
            % INPUT: none
            % OUTPUT: DROP [1x1 struct]
            %
            % drop.H = [Nx1 double]; % mean curvatures (1/vx)
            % drop.K = [Nx1 double]; % gaussian curvatures (1/vx)
            % drop.x = [Nx1 double]; % x coordinates (vx)
            % drop.y = [Nx1 double]; % y coordinates (vx)
            % drop.z = [Nx1 double]; % z coordinates (vx)
            % drop.azi = [Nx1 double]; % azimuthal coordinate (radians)
            % drop.elv = [Nx1 double]; % elevation coordinate (radians)
            % drop.r = [Nx1 double]; % radial coordinate (vx)
            % drop.lfitp = [Nx3 double]; % location fitting parameters (sigma, amplitude, offset)
            % drop.paths.tif = [1xT char]; % path to .TIF file
            % drop.paths.matfile = [1xM char]; % path to .MAT file
            %
            %   NOTES: N is the number of coordinates in the point cloud.
            %   T is the number of characters in the path to .TIF file.
            %   M is the number of characters in the path to .MAT file.
            %;
            DropRecon3Dv1.turnOffWarnings();
            DropRecon3Dv1.mexCheck();
            defPar = DropRecon3Dv1.getDefPar();
            N = defPar.numOfExpTifs;
            for n = 1:N
                if license('test', 'Distrib_Computing_Toolbox')
                    poolObj = parpool();
                end
                fprintf('\nStarting on %d of %d\n',n,N);
                defPar.tifNum = n;
                fprintf('\n Loading tiff stack. This may take a few seconds... \n')
                
                image = DropRecon3Dv1.defPar2Image(defPar);
                
                fprintf('\n Constructing point cloud. This may take a minute... \n')
                tic
                recon = DropRecon3Dv1.image2Recon(image);
                toc
                fprintf('\n Measuring surface curvatures. This may take a minute... \n')
                tic
                meas = DropRecon3Dv1.recon2Measurement(recon);
                toc
                fprintf('\n Saving analysis and shutting down parallel pool... \n')
                
                drop = DropRecon3Dv1.saveDrop(meas,defPar);
                
                if license('test', 'Distrib_Computing_Toolbox')
                    delete(poolObj)
                end
                fprintf('\n ...Finished %d of %d\n\n',n,N);
                
            end
        end
        function [grph] = drop2graphics(drop)
            close all;
            clc;
            
            n = 2e2;
            
            azi = [drop.azi-pi;drop.azi+pi];
            elv = [drop.elv;drop.elv];
            r_um = [drop.r.um;drop.r.um];
            r_um_int = scatteredInterpolant(azi,elv,r_um,'natural');
            H_inv_um = [drop.H.inv_um;drop.H.inv_um];
            H_inv_um_int = scatteredInterpolant(azi,elv,H_inv_um,'natural');
            
            agv = linspace(-pi,pi,n);
            
            egv = linspace(-pi/2,pi/2,n);
            [A,E] = meshgrid(agv,egv);
            R = r_um_int(A,E);
            H = H_inv_um_int(A,E);
            
            [X,Y,Z] = sph2cart(A,E,R);
%             
%             figure(1)
%             Ro = mean(drop.r.um);
%             r_err = (R./Ro - 1);
%             surface(X,Y,Z,r_err,'LineStyle','none')
%             axis equal
%             axis off;
%             caxis([-0.03 0.03])
%             fig = gcf;
%             fig.Color = [1 1 1];
%             
%             figure(2)
%             H_err = (H.*Ro - 1);
%             surface(X,Y,Z,H_err,'LineStyle','none')
%             axis equal
%             axis off;
%             caxis([-0.3 0.3])
%             fig = gcf;
%             fig.Color = [1 1 1];
            
            grph.X = X;
            grph.Y = Y;
            grph.Z = Z;
            grph.A = A;
            grph.E = E;
            grph.R = R;
            grph.H = H;
            
            %% MAKE RADIAL MAP GIF
            caxisStr = 'Radius (\mum)';
            gifName = 'r_map_um.gif';
            a = strfind(drop.path2ResultsFile,'/');
            b = a(end);
            gifPath = [drop.path2ResultsFile(1:b),gifName];
            DropRecon3Dv1.makeGif(X,Y,Z,R,caxisStr,gifPath)
            
            %% MAKE MEAN CURVATURE MAP GIF
            caxisStr = 'Mean Curvature (1/\mum)';
            gifName = 'H_map_inv_um.gif';
            a = strfind(drop.path2ResultsFile,'/');
            b = a(end);
            gifPath = [drop.path2ResultsFile(1:b),gifName];
            DropRecon3Dv1.makeGif(X,Y,Z,H,caxisStr,gifPath)
        end
        function [drop] = saveDrop(meas,defPar)
            drop.H.inv_px = meas.H_rec;
            drop.H.inv_um = meas.H_rec./meas.vsx;
            drop.K.inv_px = meas.K_rec;
            drop.K.inv_um = meas.K_rec./meas.vsx;
            drop.azi = meas.azi;
            drop.elv = meas.elv;
            drop.r.px = meas.r;
            drop.r.um = meas.r.*meas.vsx;
            drop.r.px = meas.r;
            drop.r.um = meas.r.*meas.vsx;
            drop.x.px = meas.x;
            drop.x.um = meas.x.*meas.vsx;
            drop.y.px = meas.y;
            drop.y.um = meas.y.*meas.vsx;
            drop.z.px = meas.z;
            drop.z.um = meas.z.*meas.vsx;
            drop.vsx = meas.vsx;
            drop.vsy = meas.vsy;
            drop.vsz = meas.vsz;
            drop.path2ResultsFile = meas.path2ResultsFile;
            drop.path2Tif = meas.path2Tif;
            traceData.widths.px = meas.fitWidth;
            drop.traceData = traceData;
            
            if ~defPar.noGifs
                 grph = DropRecon3Dv1.drop2graphics(drop);
            end
            
            save(drop.path2ResultsFile,'drop','defPar','grph');
        end
        function [defPar] = getDefPar()
            % [DEFPAR] = GETDEFPAR() Prompts user to define various
            % paramters needed for the reconstruction and analysis of the
            % drop.
            %
            % OUTPUT: defPar [1x1 struct] containing multiple fields
            % associated with the reconstruction, surface analysis, and the
            % tif(s) being analyzed.
            %
            %% Prompt user to select directory containing tifs.
            prmptStr = ['Please the select the directory containing the',...
                ' tif(s) you would like to analyze.\n'];
            fprintf(['\n',prmptStr]);
            d = uigetdir(prmptStr);
            directory = [d,'/'];
            filenameStruct = dir([d,'/*.tif']);
            numOfExpTifs = length(filenameStruct);
            fprintf('Analyzing the following tifs:')
            filenames = cell(numOfExpTifs,1);
            for n = 1:numOfExpTifs
                tifname = filenameStruct(n).name;
                filenames{n} = tifname;
                path2tif = strrep(tifname,'\','/');
                fprintf([path2tif,'\n'])
            end
            fprintf('\n');
            
            %% Prompt user for voxel sizes
            getVoxelSizeFromFile = 0;
            vsx = input('Please specify voxel size in x and y (microns): ');
            vsz = input('Please specify voxel size in z (microns): ');
            
            
            zxRatio = vsz/vsx;
            if zxRatio>3 || zxRatio<(1/3)
                warning(['It is recommended that the voxel sizes be ',...
                    'approximately equal in x,y and z\n'])
            end
            
            %% Give user the option to supply custom parameters.
            fprintf('\nPlease make a selection in the dialog box.\n')
            useDefault = ...
                questdlg('Would you like to use default parameters?',...
                'Default or Custom','Custom','Default','Default');
            
            switch useDefault
                case 'Default'
                    d_lambda = 1;
                    numOfIterations = 2;
                    nPointsPerPixel = 2;
                    maxPointsFitted = 1e5;
                    med_fltr_sz = 2;
                    strbl = 0;
                    gauss = 0;
                    useStrblFilter = 0;
                    useGaussFilter = 0;
                    findCurvs = 1;
                    noGifs = 0;
                    
                case 'Custom'
                    while true
                        switch input('Would you like to include curvature analysis? (Y/N)\n','s')
                            case {'n','N','no','No'}
                                findCurvs = 0;
                                break
                            case {'y','Y','yes','Yes'}
                                findCurvs = 1;
                                break
                            otherwise
                                fprintf('Invalid response. Try again.\n')
                        end
                    end
                    d_lambda = input(['Specify the density of points ',...
                        'at surface\n(default is 1):  ']);
                    nPointsPerPixel = input(['Specify the ',...
                        'density of resampled intensties during ray ',...
                        'tracing\n(default is 2):  ']);
                    if findCurvs
                        numOfIterations = input(['Specify the ',...
                            'number of iterations of surface fitting \n',...
                            '(default is 2):  ']);
                    else
                        numOfIterations = 2;
                    end
                    maxPointsFitted = input(['Specify the ',...
                        'number of points to fit \n',...
                        '(default is 1e5):  ']);
                    med_fltr_sz = input(['Specify median filter size\n',...
                        '(default is 2):  ']);
                    while true
                        switch input(['Enter ''G'' to use a gaussian filter,',...
                                '''S'' to use a steerable filter, or ''N''',...
                                ' to use neither (default)\n'],'s')
                            case {'G','g'}
                                useGaussFilter = 1;
                                useStrblFilter = 0;
                                strbl = 0;
                                gauss = input('Enter size of gaussian filter\n (order 1 recommended):');
                                break
                            case {'S','s'}
                                useGaussFilter = 0;
                                useStrblFilter = 1;
                                strbl = input('Enter size of steerable filter\n (order 1 recommended):');
                                gauss = 0;
                                break
                            case {'N','n'}
                                useGaussFilter = 0;
                                useStrblFilter = 0;
                                strbl = 0;
                                gauss = 0;
                                break
                            otherwise
                                fprintf('Invalid response. Try again.\n')
                        end
                    end
                    while true
                        switch input('Would you like to make animated gifs? (Y/N)\n','s')
                            case {'n','N','no','No'}
                                noGifs = 1;
                                break
                            case {'y','Y','yes','Yes'}
                                noGifs = 0;
                                break
                            otherwise
                                fprintf('Invalid response. Try again.\n')
                        end
                    end
            end
            
            %% Other Parameters
            IQRcoordFilter = 0;
            patchMethod = 'onePxNormChng';
            beta = 2;
            degPoly = 2;
            N = beta*3^degPoly;
            minRp = sqrt(N/pi)*d_lambda;
            maxHVar = 1;
            shuffleRays = 0;
            ellipticalPatch = 0;
            chiSqrdTest = 1;
            tifNum = 1;
            numOfModelTifs = 0;
            
            %% Assign parameters as fields of defPar
            
            % DIRECTORY
            defPar.directory = directory;
            defPar.numOfExpTifs = numOfExpTifs;
            defPar.filenames = filenames;
            defPar.useDefault = useDefault;
            
            % VOXEL SIZES
            defPar.getVoxelSizeFromFile = getVoxelSizeFromFile;
            defPar.vsx = vsx;
            defPar.vsy = defPar.vsx;
            defPar.vsz = vsz;
            defPar.zxRatio = zxRatio;
            
            % CUSTOMIZABLE PARAMETERS
            defPar.d_lambda = d_lambda;
            defPar.numOfIterations = numOfIterations;
            defPar.nPointsPerPixel = nPointsPerPixel;
            defPar.maxPointsFitted = maxPointsFitted;
            defPar.medianFilterSizeInPixels = med_fltr_sz;
            filterSize.strbl = strbl;
            filterSize.gauss = gauss;
            defPar.filterSize = filterSize;
            defPar.noGifs = noGifs; % 0 means no animated gifs
            defPar.IQRcoordFilter = IQRcoordFilter;
            defPar.useGaussFilter = useGaussFilter;
            defPar.useStrblFilter = useStrblFilter;
            defPar.findCurvs = findCurvs;
            
            % OTHER PARAMETERS
            defPar.patchMethod = patchMethod;
            defPar.beta = beta;
            defPar.degPoly = degPoly;
            defPar.minRp = minRp;
            defPar.maxHVar = maxHVar;
            defPar.shuffleRays = shuffleRays;
            defPar.ellipticalPatch = ellipticalPatch;
            defPar.chiSqrdTest = chiSqrdTest;
            defPar.tifNum = tifNum;
            defPar.medFiltPatch.idx = [];
            defPar.H_rec = [];
            defPar.numOfModelTifs = numOfModelTifs;
            
        end
        function [image] = defPar2Image(defPar)
            % [IMAGE] = DEFPAR2IMAGE(DEFPAR) Uses parameters defined by
            % DEFPAR to load and process image before reconstruction and
            % analysis of point cloud.
            %
            % OUTPUT: image [1x1 struct] containing fields associated with
            % the image being loaded into the workspace and processed for
            % segmentation.
            %;
            image = defPar;
            image.filename = defPar.filenames{defPar.tifNum};
            image.path2Tif = [defPar.directory,image.filename];
            image.rawInt = DropRecon3Dv1.LoadTif(image.path2Tif);
            
            image.path2ResultsDir = [defPar.directory,...
                image.filename(1:end-4),'-Results/'];
            if ~exist(image.path2ResultsDir,'dir')
                mkdir(image.path2ResultsDir);
            end
            image.path2ResultsFile = [image.path2ResultsDir,...
                defPar.filenames{defPar.tifNum}(1:end-4),'.mat'];
            
            [image.nPixelsX,image.nPixelsY,image.nPixelsZ]=...
                size(image.rawInt); % determines size of the tif stack
            % parameters to define aspect ratio and center of ellipsoid in pixel units.
            
            [X,Y,Z] = ndgrid(1:1:image.nPixelsX,1:1:image.nPixelsY,1:1:image.nPixelsZ); % creates mesh
            % of size of filtered tif stack
            
            if image.medianFilterSizeInPixels>0
                image = DropRecon3Dv1.applyMedianFilter2Image(image);
                V=double(image.medFltInt); % assigns filtered tif stack volume to 'V'
            else
                V=double(image.rawInt); % assigns filtered tif stack volume to 'V'
            end
            image.getRawInt = griddedInterpolant(X,Y,Z,V,'nearest');
            if image.useStrblFilter
                M = 2;% Surface Detector
                [V_filtered,~,~]=steerableDetector3D(V,M,...
                    image.filterSize.strbl,image.zxRatio);
                image.fInterp = griddedInterpolant(X,Y,Z,V_filtered,'cubic');%,'none'); % Constructs new slices of
                % "drop" in the volume by interpolation so that mesh is uniformly filled.
            else
                image.fInterp = griddedInterpolant(X,Y,Z,V,'cubic');%,'none'); % Constructs new slices of
                % "drop" in the volume by interpolation so that mesh is uniformly filled.
            end
            [YI,XI,ZI] = meshgrid(1:1:image.nPixelsY,1:1:image.nPixelsX,1:1/image.zxRatio:image.nPixelsZ); %
            image.interpDbl=image.fInterp(XI,YI,ZI);
            
            %[~,~,nVirtualStacks] = size(ZI);
            
            if image.useGaussFilter
                V = image.interpDbl;
                image.gaussPatchSize = 2*(floor((9*image.filterSize.gauss)/2)+0.5);
                image.interpDbl = smooth3(V,'gaussian',...
                    image.gaussPatchSize,image.filterSize.gauss);
            end
            image.interp16Bit = uint16(image.interpDbl);
            
        end
        function [intArray3D] = LoadTif(pathtotif)
            FileTif=pathtotif;
            InfoImage=imfinfo(FileTif);
            mImage=InfoImage(1).Width;
            nImage=InfoImage(1).Height;
            NumberImages=length(InfoImage);
            intArray3D=zeros(nImage,mImage,NumberImages,'uint16');
            TifLink = Tiff(FileTif, 'r');
            for i=1:NumberImages
                TifLink.setDirectory(i);
                intArray3D(:,:,i)=TifLink.read();
            end
            TifLink.close();
        end
        function [image] = applyMedianFilter2Image(image)
            if image.medianFilterSizeInPixels
                [rows,cols,pages] = size(image.rawInt);
                imIn = image.rawInt;
                imOut = uint16(zeros(rows,cols,pages));
                Lxy = round(image.medianFilterSizeInPixels);
                Lz = round(Lxy.*image.vsx);
                for i = 1:rows
                    rowMin = i-Lxy;
                    if rowMin<1
                        rowMin=1;
                    end
                    rowMax = i+Lxy;
                    if rowMax>rows
                        rowMax=rows;
                    end
                    for j=1:cols
                        colMin = j-Lxy;
                        if colMin<1
                            colMin=1;
                        end
                        colMax = j+Lxy;
                        if colMax>cols
                            colMax=cols;
                        end
                        for k=1:pages
                            pagMin = k-Lz;
                            if pagMin<1
                                pagMin=1;
                            end
                            pagMax = k+Lz;
                            if pagMax>pages
                                pagMax=pages;
                            end
                            imBox = imIn(rowMin:rowMax,colMin:colMax,pagMin:pagMax);
                            imOut(i,j,k) = median(imBox(:));
                        end
                    end
                end
                image.medFltInt = imOut;
            else
                image.medFltInt = image.rawInt;
            end
        end
        function [recon] = image2Recon(image)
            %;
            recon = image;
            recon.oversamplingFactor = 1;
            recon.fitAlongNormals = 1;
            recon = DropRecon3Dv1.toss(recon);
            recon = DropRecon3Dv1.reconCoords(recon);
            [recon] = DropRecon3Dv1.applyMedianFilter2Recon(recon);
            [recon] = DropRecon3Dv1.findAxiallySymmetricEllipsoid(recon);
            [recon] = DropRecon3Dv1.findGenEllipsoid(recon);
            [recon] = DropRecon3Dv1.makeGriddedSurface(recon);
            
            
        end
        function [recon] = toss(recon)
            
            FitAlongNormals = 0; % just for the quick one
            SingleDrop = recon.interp16Bit;
            d_lambda = recon.d_lambda;
            NPointsPerPixel = recon.nPointsPerPixel;
            acceptRange = recon.IQRcoordFilter;
            [NPixelsX, NPixelsY, NPixelsZ]=size(SingleDrop);
            
            %1.A interpolate image in 3D linearly
            [Y,X,Z] = meshgrid(1:1:NPixelsY,1:1:NPixelsX,1:1:NPixelsZ);
            V=double(SingleDrop);
            F = griddedInterpolant(X,Y,Z,V,'linear','none');
            G = griddedInterpolant(X,Y,Z,V,'nearest','none');
            %%  first guess of center-of-mass and radius
            
            %2.A Center of mass of the interface
            recon.threshold1A=0.7;
            COM=DropRecon3Dv1.getCenterOfMass(SingleDrop,recon.threshold1A);
            
            %2.B center of mass of drop by maximum amplitude projection
            recon.threshold1B=0.05;
            [centerY1, centerZ1,radiusY1,radiusZ1]=...
                DropRecon3Dv1.getCenterOfMassMIP(SingleDrop,1,round(COM),...
                recon.threshold1B);
            [centerX1, centerZ2,radiusX1,radiusZ2]=...
                DropRecon3Dv1.getCenterOfMassMIP(SingleDrop,2,round(COM),...
                recon.threshold1B);
            [centerX2, centerY2,radiusX2,radiusY2]=...
                DropRecon3Dv1.getCenterOfMassMIP(SingleDrop,3,round(COM),...
                recon.threshold1B);
            
            % average results
            centerX=(centerX1+centerX2)/2;
            centerY=(centerY1+centerY2)/2;
            centerZ=(centerZ1+centerZ2)/2;
            radiusX=(radiusX1+radiusX2)/2;
            radiusY=(radiusY1+radiusY2)/2;
            radiusZ=(radiusZ1+radiusZ2)/2;
            radii=[radiusX,radiusY,radiusZ];
            
            %%   coarse determination of interface coordinates
            
            %    generate coarse mesh of sampling points on a sphere
            CoarseRayNumber = 100;
            H_dummy = ones(CoarseRayNumber,1);
            [RV1,~,~]=DropRecon3Dv1.getSpiraledRaysEllipsoidMeanNnDist(...
                d_lambda,max(radii),max(radii),max(radii),CoarseRayNumber,recon.shuffleRays);
            d_lambdaCoarse = sqrt((4*pi*mean(radii)^2)/CoarseRayNumber);
            
            NPoints=length(RV1);
            RStart1=[centerX*ones(NPoints,1),centerY*ones(NPoints,1),...
                centerZ*ones(NPoints,1)];
            fltrPatchSize = recon.medianFilterSizeInPixels;
            [CoarseCoords,~,~,~] = ...
                DropRecon3Dv1.reconstructCoords(RStart1,RV1,F,FitAlongNormals,...
                H_dummy,NPointsPerPixel,acceptRange,fltrPatchSize);
            
            centerX = mean(CoarseCoords(:,1));
            centerY = mean(CoarseCoords(:,2));
            centerZ = mean(CoarseCoords(:,3));
            COM = [centerX,centerY,centerZ];
            
            [numOfCoords,~] = size(CoarseCoords);
            CenteredCoarseCoords = CoarseCoords-ones(numOfCoords,1)*COM;
            if numOfCoords ~= CoarseRayNumber
                warning('Not all rays "thrown" were "caught" at interface...')
            end
            
            [~,eigenVectors,~] = DropRecon3Dv1.getEigen(CenteredCoarseCoords);
            
            principalAxis_1 = eigenVectors(:,1);
            principalAxis_2 = eigenVectors(:,2);
            principalAxis_3 = eigenVectors(:,3);
            
            pa1_x = principalAxis_1(1);
            pa1_y = principalAxis_1(2);
            pa1_z = principalAxis_1(3);
            pa1_r = sqrt(pa1_x^2+pa1_y^2);
            refPolarAngle = atan(abs(pa1_r/pa1_z));
            PolarAngle = refPolarAngle*(pa1_z>=0)+...
                (pi-refPolarAngle)*(pa1_z<0);
            refAziAngle = atan(abs(pa1_y/pa1_x));
            AzimuthalAngle = refAziAngle*(pa1_x>=0)*(pa1_y>=0)+...
                (pi-refAziAngle)*(pa1_x<0)*(pa1_y>=0)+...
                (pi+refAziAngle)*(pa1_x<0)*(pa1_y<0)+...
                (2*pi-refAziAngle)*(pa1_x>=0)*(pa1_y<0);
            
            CenteredCoords_rot = DropRecon3Dv1.rotateCoordinates(...
                CenteredCoarseCoords,-PolarAngle,-AzimuthalAngle);
            
            CenteredCoordsX_rot = CenteredCoords_rot(:,1);
            CenteredCoordsY_rot = CenteredCoords_rot(:,2);
            CenteredCoordsZ_rot = CenteredCoords_rot(:,3);
            
            radiusX1 = max(CenteredCoordsX_rot);
            radiusX2 = abs(min(CenteredCoordsX_rot));
            
            radiusY1 = max(CenteredCoordsY_rot);
            radiusY2 = abs(min(CenteredCoordsY_rot));
            
            radiusZ1 = max(CenteredCoordsZ_rot);
            radiusZ2 = abs(min(CenteredCoordsZ_rot));
            
            radiusX = (radiusX1+radiusX2)/2;
            radiusY = (radiusY1+radiusY2)/2;
            radiusZ = (radiusZ1+radiusZ2)/2;
            
            radii=[radiusX,radiusY,radiusZ];
            
            %%   generate fine mesh of sampling points on a sphere
            [RV1,meanNNDist,samplingDens,H_est,K_est]=...
                DropRecon3Dv1.getSpiraledRaysEllipsoidMeanNnDist(...
                d_lambda,radii(1),radii(2),radii(3));
            RV1 = DropRecon3Dv1.rotateCoordinates(RV1,PolarAngle,AzimuthalAngle);
            
            recInfo = {};
            recInfo.meanNNDist = meanNNDist;
            recInfo.samplingDens = samplingDens;
            
            NPoints=length(RV1);
            RStart1=[centerX*ones(NPoints,1),centerY*ones(NPoints,1),...
                centerZ*ones(NPoints,1)];
            
            % Pass values to structured variable 'recon'
            recon.rayStart = RStart1;
            recon.traceVectors = RV1;
            recon.fInterp = F;
            recon.recInfo = recInfo;
            recon.roughHEst = H_est;
            recon.roughKEst = K_est;
            recon.gInterp = G;
        end
        function [COM] = getCenterOfMass(image3D,threshold)
            
            [NPixelsX,~,NPixelsZ]=size(image3D);
            maxI=max(image3D(:));
            ts=threshold*maxI;
            
            %threshold the image
            image3Dts=image3D.*uint16(image3D>ts);
            
            Coords=[];
            Int=[];
            for k=1:NPixelsZ
                singlePlane=squeeze(image3Dts(:,:,k));
                indInt=find(singlePlane>0);
                [X,Y] = ind2sub(NPixelsX,indInt) ;
                Coords=[Coords; [X Y k*ones(length(indInt),1)]];
                Int=[Int; singlePlane(indInt)];
            end
            weightedCoordsX=Coords(:,1).*double(Int);
            weightedCoordsY=Coords(:,2).*double(Int);
            weightedCoordsZ=Coords(:,3).*double(Int);
            
            IntSum=sum(Int);
            %NPixelsLowR=NPixelsXLowR*NPixelsYLowR*NPixelsZLowR;
            COM=[sum(weightedCoordsX)/IntSum,sum(weightedCoordsY)/IntSum,sum(weightedCoordsZ)/IntSum];
        end
        function [Center1,Center2,Radius1,Radius2]=getCenterOfMassMIP(...
                image3D,direction,initialCenterCoords,threshold )
            if isnan(sum(initialCenterCoords))
                if sum(image3D(:))==0
                    error('Image array is all zeros.')
                else
                    error('Problem with center coordinates passed to DropRecon3Dv1.getCenterOfMassMIP.')
                end
            end
            initialCenterCoordX=initialCenterCoords(1);
            initialCenterCoordY=initialCenterCoords(2);
            initialCenterCoordZ=initialCenterCoords(3);
            
            [NPixelsX,NPixelsY,NPixelsZ]=size(image3D);
            
            maxIp=max(image3D,[],direction);
            
            if direction==1
                maxIp=reshape(maxIp,NPixelsY,NPixelsZ);
                maxAmp=max(max(maxIp));
                
                %find center in Y
                
                %start at initial guessed center, increases Y, looks for pixel at wich Intensity drops below threshold
                %indVolume: Offset in pixel with respect to initial guess center
                find(maxIp(initialCenterCoordY:1:NPixelsY,initialCenterCoordZ)>threshold*maxAmp,1,'last');
                indVolumePlus1=find(maxIp(initialCenterCoordY:1:NPixelsY,initialCenterCoordZ)>threshold*maxAmp,1,'last');
                
                %start at initial guessed center, decreases Y, looks for pixel at wich Intensity drops below threshold
                indVolumeMinus1=find(maxIp(initialCenterCoordY:-1:1,initialCenterCoordZ)>threshold*maxAmp,1,'last');
                Center1=(indVolumePlus1-indVolumeMinus1)/2+initialCenterCoordY;
                Radius1=(indVolumePlus1+indVolumeMinus1)/2;
                
                %find center in Z
                indVolumePlus2=find(maxIp(initialCenterCoordY,initialCenterCoordZ:1:NPixelsZ)>threshold*maxAmp,1,'last');
                indVolumeMinus2=find(maxIp(initialCenterCoordY,initialCenterCoordZ:-1:1)>threshold*maxAmp,1,'last');
                Center2=(indVolumePlus2-indVolumeMinus2)/2+initialCenterCoordZ;
                Radius2=(indVolumePlus2+indVolumeMinus2)/2;
                
            end;
            
            if direction==2
                
                maxIp=reshape(maxIp,NPixelsX,NPixelsZ);
                maxAmp=max(max(maxIp));
                
                %find center in X
                indVolumePlus1=find(maxIp(initialCenterCoordX:NPixelsX,initialCenterCoordZ)>threshold*maxAmp,1,'last');
                indVolumeMinus1=find(maxIp(initialCenterCoordX:-1:1,initialCenterCoordZ)>threshold*maxAmp,1,'last');
                Center1=(indVolumePlus1-indVolumeMinus1)/2+initialCenterCoordX;
                Radius1=(indVolumePlus1+indVolumeMinus1)/2;
                %find center in Z
                indVolumePlus2=find(maxIp(initialCenterCoordX,initialCenterCoordZ:NPixelsZ)>threshold*maxAmp,1,'last');
                indVolumeMinus2=find(maxIp(initialCenterCoordX,initialCenterCoordZ:-1:1)>threshold*maxAmp,1,'last');
                Center2=(indVolumePlus2-indVolumeMinus2)/2+initialCenterCoordZ;
                Radius2=(indVolumePlus2+indVolumeMinus2)/2;
            end;
            
            if direction==3
                maxAmp=max(max(maxIp));
                maxIp=reshape(maxIp,NPixelsX,NPixelsY);
                
                %find center in X
                indVolumePlus1=find(maxIp(initialCenterCoordX:NPixelsX,initialCenterCoordY)>threshold*maxAmp,1,'last');
                indVolumeMinus1=find(maxIp(initialCenterCoordX:-1:1,initialCenterCoordY)>threshold*maxAmp,1,'last');
                Center1=(indVolumePlus1-indVolumeMinus1)/2+initialCenterCoordX;
                Radius1=(indVolumePlus1+indVolumeMinus1)/2;
                %find center in Y
                indVolumePlus2=find(maxIp(initialCenterCoordX,initialCenterCoordY:NPixelsY)>threshold*maxAmp,1,'last');
                indVolumeMinus2=find(maxIp(initialCenterCoordX,initialCenterCoordY:-1:1)>threshold*maxAmp,1,'last');
                Center2=(indVolumePlus2-indVolumeMinus2)/2+initialCenterCoordY;
                Radius2=(indVolumePlus2+indVolumeMinus2)/2;
                
            end
            
        end
        function [Rays,meanNNDist,samplingDens,H_est,K_est] = ...
                getSpiraledRaysEllipsoidMeanNnDist(d_lambda,a,b,c,...
                numberOfRays,ShuffleRays)
            if nargin~=6
                ShuffleRays=0;% 1 => Do shuffle rays, 0 => Do not shuffle rays
            end
            if ShuffleRays==0
                %fprintf('\nRays will NOT be shuffled...\n')
            else
                %fprintf('\nRays will be shuffled...\n')
            end
            %fprintf('\n Runnning getSpiraledRaysEllipsoidMeanNnDist...\n')
            
            surfaceArea = DropRecon3Dv1.getSurfaceAreaEllipsoid( a,b,c );
            if nargin<5
                numberOfRays = round(surfaceArea./(d_lambda.^2));
            end
            
            z = linspace(1-1/numberOfRays,1/numberOfRays-1,numberOfRays);
            radius=sqrt(1-z.^2);
            
            goldenAngle = pi*(3-sqrt(5));
            theta = goldenAngle*(1:numberOfRays);
            
            Rays = zeros(numberOfRays,3);
            Rays(:,1) = a*radius.*cos(theta);
            Rays(:,2) = b*radius.*sin(theta);
            Rays(:,3) = c*z;
            %
            %Rays(:,1) = radius.*cos(theta);
            %Rays(:,2) = radius.*sin(theta);
            %Rays(:,3) = z;
            %polAngle = pi*rand(1);
            %aziAngle = 2*pi*rand(1);
            %Rays = DropRecon3Dv1.rotateCoordinates(Rays, polAngle, aziAngle);
            
            nnd = DropRecon3Dv1.getNNDistances(Rays);
            meanNNDist = mean(nnd);
            samplingDens = numberOfRays/surfaceArea;
            
            % Shuffle Rays;
            if ShuffleRays
                ind = randperm(length(Rays));
                Rays = Rays(ind,:);
            end
            %% ESTIMATE CURVATURE FROM ECCENTRICITIES
            % u is [0, 2pi], v is [0, pi]
            x = Rays(:,1);
            y = Rays(:,2);
            z = Rays(:,3);
            
            cosv = z./c;
            v = acos(cosv);
            sinv = sin(v);
            cosu = x./(a*sinv);
            u = acos(cosu).*(y>=0)+(2*pi - acos(cosu)).*(y<0);
            sinu = sin(u);
            
            H_est = a*b*c*(3*(a^2+b^2)+2*c^2+(a^2+b^2-2*c^2)*cos(2*v)-...
                2*(a^2-b^2)*cos(2*u).*sinv.^2)./...
                (8*(a^2*b^2*cosv.^2+c^2*(b^2*cosu.^2+...
                a^2*sinu.^2).*sinv.^2).^(3/2));
            K_est = ((a*b*c)^2)./((a*b*cosv).^2+(c^2).*((b*cosu).^2+...
                (a*sinu).^2).*sinv.^2).^2;
            
            %fprintf('\n ...Done!\n')
        end
        function [surfaceArea] = getSurfaceAreaEllipsoid(a,b,c)
            %http://www.web-formulas.com/Math_Formulas/Geometry_Surface_of_Ellipsoid.aspx
            p=1.6075;
            surfaceArea=4*pi*((a^p*b^p+a^p*c^p+b^p*c^p)/3)^(1/p);
        end
        function [NNDistances,NNInd]=getNNDistances(CoordsIn,N)
            NS_kd = KDTreeSearcher(CoordsIn);
            if nargin < 2
                [ind, nd] = knnsearch(NS_kd,CoordsIn,'k',2);
                clear NS_kd;
                [m,n] = size(ind);
                if n<2 || m<1
                    NNDistances = NaN;
                    NNInd = NaN;
                else
                    NNDistances=nd(:,2);
                    NNInd=ind(:,2);
                end
            else
                [ind, nd] = knnsearch(NS_kd,CoordsIn,'k',N+1);
                NNDistances=mean(nd(:,2:N+1),2);
                NNInd=ind(:,2:N+1);
            end
        end
        function [CoordsRec,ErrCoordsRec,InterfaceWidthFromFit,H_out,...
                RayTraceInfo] = reconstructCoords(RStart1,RV1,F,...
                FitAlongNormals,H_est,NPointsPerPixel,acceptRange,...
                FltrPatchSize)
            
            
            
            %% 2.C fit Gaussians to determine interface
            
            %generateScatterPlot(RV1(:,1),RV1(:,2),RV1(:,3),'trace rays')
            
            %parameters for fit algorithm, do not change
            threshold1=0.4;threshold2=0.05;
            
            %hard coded, works for most drops
            rVscale1=6;    %   initial ray is 6x larger than guessed radius
            rVscale2=0.5;   %	refined ray is centered around and 0.5 x length of guessed radius
            
            [fCoords,errfCoords,RayTraceInfo,~,sd_est,H_out]=...
                DropRecon3Dv1.getFittedCoordsCPP(RStart1,RV1,NPointsPerPixel,...
                F,threshold1,threshold2,rVscale1,rVscale2,H_est);
            
            %generateScatterPlot(fCoords(:,1),fCoords(:,2),fCoords(:,3),'coords first recontr step')
            %figure();plot(errfCoords)
            
            %   final determination of COM and radius
            
            COMFit=mean(fCoords);
            %fCoordsCOM=(fCoords-ones(length(fCoords),1)*COMFit);
            %RadiusFit=mean(sqrt(fCoordsCOM(:,1).^2+fCoordsCOM(:,2).^2+fCoordsCOM(:,3).^2));
            %DistFromCOM=sqrt(fCoordsCOM(:,1).^2+fCoordsCOM(:,2).^2+fCoordsCOM(:,3).^2);
            %hist(DistFromCOM,1000)
            
            %%   optimize sampling
            
            % create more rays and
            
            %%  fit interface coords (radially traced)
            fCoordsR=fCoords;
            
            % remove NaN coordinates
            indFiltered = sum(isnan(fCoords),2)==0;
            fCoordsR = fCoordsR(indFiltered,:);
            errfCoordsR.width = errfCoords.width(indFiltered);
            errfCoordsR.center = errfCoords.center(indFiltered);
            errfCoordsR.amp = errfCoords.amp(indFiltered);
            errfCoordsR.offset = errfCoords.offset(indFiltered);
            
            H_out = H_out(indFiltered);
            RV1 = RV1(indFiltered,:);
            sd_est = sd_est(indFiltered);
            if isempty(fCoordsR)
                CoordsRec = [];
                ErrCoordsRec = [];
                H_out = [];
                InterfaceWidthNormalFit = [];
                RayTraceInfo.rawIntsAtPeak = [];
                RayTraceInfo.fitIntsAtPeak = [];
                RayTraceInfo.fitOffset = [];
            else
                
                if acceptRange>0
                    %NN_IQF = 1;
                    ECR_Width_IQF = 1;
                    ECR_Center_IQF = 0;
                    ECR_Amp_IQF = 0;
                    ECR_Offset_IQF = 0;
                else
                    %NN_IQF = 0;
                    ECR_Width_IQF = 0;
                    ECR_Center_IQF = 0;
                    ECR_Amp_IQF = 0;
                    ECR_Offset_IQF = 0;
                end
                % OLD - FIRST use InterQuartile Filter to remove major outliers
                NN_MDF = 1;
                %% FIRST use a Minimum Density Filter to throw out points
                % that are isolated by considering the number in the neighborhood
                if NN_MDF
                    
                    patchSize = FltrPatchSize;
                    numInNeighborhood = zeros(length(fCoordsR),1);
                    NS_kd = KDTreeSearcher(fCoordsR);
                    if length(fCoordsR)>4.0e4 && license('test', 'Distrib_Computing_Toolbox') % 19 Aug 2016
                        
                        parfor i = 1:length(fCoordsR)
                            %% GET CIRCULAR PATCH INDICES
                            queryPoint = fCoordsR(i,:);
                            [idx,~]=rangesearch(NS_kd,queryPoint,patchSize); % this delivers the indices of points in the neighborhood
                            idx = idx{1};
                            idx=nonzeros(idx);
                            numInNeighborhood(i) = length(idx);
                        end
                        
                    else
                        for i = 1:length(fCoordsR)
                            %% GET CIRCULAR PATCH INDICES
                            queryPoint = fCoordsR(i,:);
                            [idx,~]=rangesearch(NS_kd,queryPoint,patchSize); % this delivers the indices of points in the neighborhood
                            idx = idx{1};
                            idx=nonzeros(idx);
                            %this is the neighborhood
                            %patchCoords=fCoordsR(idx,:);
                            %[numInNeighborhood(i),~] = size(patchCoords);
                            numInNeighborhood(i) = length(idx);
                        end
                    end
                    
                    %indFiltered_NN = numInNeighborhood>(0.5*pi*(patchSize/d_lambda).^2);
                    indFiltered_NN = numInNeighborhood>=(mean(numInNeighborhood)-3*std(numInNeighborhood));
                    
                else
                    indFiltered_NN = ones(length(fCoordsR),1);
                end
                %% SECOND use InterQuartile Filter to remove major outliers by
                % by the error in the width parameter of the gaussian fit
                if ECR_Width_IQF
                    indFiltered_ErrWidth = DropRecon3Dv1.filterInterquartileRange(...
                        errfCoordsR.width,acceptRange);
                else
                    indFiltered_ErrWidth = ones(length(fCoordsR),1);
                end
                %% THIRD use InterQuartile Filter to remove major outliers
                % by the error in the center parameter of the gaussian fit
                if ECR_Center_IQF
                    indFiltered_ErrCenter = DropRecon3Dv1.filterInterquartileRange(...
                        errfCoordsR.center,acceptRange);
                else
                    indFiltered_ErrCenter = ones(length(fCoordsR),1);
                end
                
                %% FOURTH use InterQuartile Filter to remove major outliers
                % by the error in the amplitude parameter of the gaussian fit
                if ECR_Amp_IQF
                    indFiltered_ErrAmp = DropRecon3Dv1.filterInterquartileRange(...
                        errfCoordsR.amp,acceptRange);
                else
                    indFiltered_ErrAmp = ones(length(fCoordsR),1);
                end
                
                %% FIFTH use InterQuartile Filter to remove major outliers
                % by error in the offset parameter of the gaussian fit.
                if ECR_Offset_IQF
                    indFiltered_ErrOffset = DropRecon3Dv1.filterInterquartileRange(...
                        errfCoordsR.offset,acceptRange);
                else
                    indFiltered_ErrOffset = ones(length(fCoordsR),1);
                end
                
                %size(indFiltered_NN);
                %size(indFiltered_ErrWidth);
                %size(indFiltered_ErrCenter);
                %size(indFiltered_ErrAmp);
                %size(indFiltered_ErrOffset);
                
                indFiltered = (indFiltered_NN&indFiltered_ErrWidth&...
                    indFiltered_ErrCenter&...
                    indFiltered_ErrAmp&...
                    indFiltered_ErrOffset);
                
                CoordsRec = fCoordsR(indFiltered,:);
                H_out = H_out(indFiltered);
                ErrCoordsRec.width = errfCoordsR.width(indFiltered);
                ErrCoordsRec.center = errfCoordsR.center(indFiltered);
                ErrCoordsRec.amp = errfCoordsR.amp(indFiltered);
                ErrCoordsRec.offset = errfCoordsR.offset(indFiltered);
                InterfaceWidthFit = sd_est(indFiltered);
                RV1 = RV1(indFiltered,:);
                RayTraceInfo.rawIntsAtPeak = RayTraceInfo.rawIntsAtPeak(indFiltered);
                RayTraceInfo.fitIntsAtPeak = RayTraceInfo.fitIntsAtPeak(indFiltered);
                RayTraceInfo.fitOffset = RayTraceInfo.fitOffset(indFiltered);
                
                if FitAlongNormals
                    %%
                    minPoints = 30;
                    %kd = KDTree(fCoordsR(indFiltered,:));
                    NS_kd = KDTreeSearcher(CoordsRec);
                    %take only the minPoints nearest neighbors, regardsless of their distance
                    [Normals,~,~] = ...
                        DropRecon3Dv1.getNormals(CoordsRec,...
                        CoordsRec,NS_kd,minPoints);%,patchRadius);
                    %orient them to point outwards
                    NormalsO = DropRecon3Dv1.setNormalOrientationRV(Normals, RV1);
                    %%  7. fit along normal direction
                    
                    nSigma=8; %length of ray across interface in units of
                    % standard deviations of the Gaussian interface
                    sd_est = median(InterfaceWidthFit);
                    rayLength=sd_est*nSigma; %length of ray across
                    % interface in units of voxel length
                    if rayLength>mean(2*sqrt(sum((fCoordsR-16).^2,2)))
                        rayLength=mean(2*sqrt(sum((fCoordsR-16).^2,2)));
                    end
                    RayDirections=NormalsO;
                    [RayStarts,Rays]=DropRecon3Dv1.getRays(rayLength,...
                        RayDirections,CoordsRec);
                    
                    %parameters for fit algorithm, do not change
                    threshold1=0.4;threshold2=0.05;
                    rVscale1=1;rVscale2=1;% scale the rays by this lengths
                    
                    [fCoordsN, errfCoordsN,RayTraceInfo,~,...
                        InterfaceWidthNormalFit,H_out]=...
                        DropRecon3Dv1.getFittedCoordsCPP(RayStarts,Rays,...
                        NPointsPerPixel,F,threshold1,threshold2,...
                        rVscale1,rVscale2,H_out);
                    
                    % FILTER OUT NAN IN fCoordsN
                    indFiltered = ~isnan(fCoordsN(:,1));
                    fCoordsN=fCoordsN(indFiltered,:);
                    errfCoordsN.width=errfCoordsN.width(indFiltered);
                    errfCoordsN.center=errfCoordsN.center(indFiltered);
                    errfCoordsN.amp=errfCoordsN.amp(indFiltered);
                    errfCoordsN.offset=errfCoordsN.offset(indFiltered);
                    H_out = H_out(indFiltered);
                    InterfaceWidthNormalFit = InterfaceWidthNormalFit(indFiltered);
                    RayTraceInfo.rawIntsAtPeak = RayTraceInfo.rawIntsAtPeak(indFiltered);
                    RayTraceInfo.fitIntsAtPeak = RayTraceInfo.fitIntsAtPeak(indFiltered);
                    RayTraceInfo.fitOffset = RayTraceInfo.fitOffset(indFiltered);
                    
                    %% Apply a minimum density filter
                    if NN_MDF
                        numOfCoords = length(fCoordsN);
                        estimatedTime = 1.2e-3*numOfCoords;
                        patchSize = FltrPatchSize;
                        numInNeighborhood = zeros(numOfCoords,1);
                        NS_kd = KDTreeSearcher(fCoordsN);
                        if estimatedTime>45 && license('test', 'Distrib_Computing_Toolbox')
                            parfor i = 1:numOfCoords
                                %% GET CIRCULAR PATCH INDICES
                                queryPoint = fCoordsN(i,:);
                                [idx,~]=rangesearch(NS_kd,queryPoint,patchSize); % this delivers the indices of points in the neighborhood
                                idx = idx{1};
                                idx=nonzeros(idx);
                                %this is the neighborhood
                                %patchCoords=fCoordsR(idx,:);
                                %[numInNeighborhood(i),~] = size(patchCoords);
                                numInNeighborhood(i) = length(idx);
                            end
                        else
                            %
                            for i = 1:numOfCoords
                                %% GET CIRCULAR PATCH INDICES
                                queryPoint = fCoordsN(i,:);
                                [idx,~]=rangesearch(NS_kd,queryPoint,patchSize); % this delivers the indices of points in the neighborhood
                                idx = idx{1};
                                idx=nonzeros(idx);
                                %this is the neighborhood
                                %patchCoords=fCoordsR(idx,:);
                                %[numInNeighborhood(i),~] = size(patchCoords);
                                numInNeighborhood(i) = length(idx);
                            end
                            %
                        end
                        %indFiltered_NN = numInNeighborhood>(0.5*pi*(patchSize/d_lambda).^2);
                        indFiltered_NN = numInNeighborhood>=(mean(numInNeighborhood)-3*std(numInNeighborhood));
                        
                    else
                        indFiltered_NN = ones(length(fCoordsN),1);
                    end
                    
                    %% SECOND use InterQuartile Filter to remove major outliers by
                    % by the error in the width parameter of the gaussian fit
                    if ECR_Width_IQF
                        indFiltered_ErrWidth = DropRecon3Dv1.filterInterquartileRange(...
                            errfCoordsN.width,acceptRange);
                    else
                        indFiltered_ErrWidth = ones(length(fCoordsN),1);
                    end
                    %% THIRD use InterQuartile Filter to remove major outliers
                    % by the error in the center parameter of the gaussian fit
                    if ECR_Center_IQF
                        indFiltered_ErrCenter = DropRecon3Dv1.filterInterquartileRange(...
                            errfCoordsN.center,acceptRange);
                    else
                        indFiltered_ErrCenter = ones(length(fCoordsN),1);
                    end
                    %% FOURTH use InterQuartile Filter to remove major outliers
                    % by the error in the amplitude parameter of the gaussian fit
                    if ECR_Amp_IQF
                        indFiltered_ErrAmp = DropRecon3Dv1.filterInterquartileRange(...
                            errfCoordsN.amp,acceptRange);
                    else
                        indFiltered_ErrAmp = ones(length(fCoordsN),1);
                    end
                    %% FIFTH use InterQuartile Filter to remove major outliers
                    % by error in the offset parameter of the gaussian fit.
                    if ECR_Offset_IQF
                        indFiltered_ErrOffset = DropRecon3Dv1.filterInterquartileRange(...
                            errfCoordsN.offset,acceptRange);
                    else
                        indFiltered_ErrOffset = ones(length(fCoordsN),1);
                    end
                    indFiltered = indFiltered_NN&indFiltered_ErrWidth&...
                        indFiltered_ErrCenter&...
                        indFiltered_ErrAmp&...
                        indFiltered_ErrOffset;
                    
                    CoordsRec=fCoordsN(indFiltered,:);
                    H_out = H_out(indFiltered);
                    ErrCoordsRec.width=errfCoordsN.width(indFiltered);
                    ErrCoordsRec.center=errfCoordsN.center(indFiltered);
                    ErrCoordsRec.amp=errfCoordsN.amp(indFiltered);
                    ErrCoordsRec.offset=errfCoordsN.offset(indFiltered);
                    InterfaceWidthFromFit = InterfaceWidthNormalFit(indFiltered);
                    RayTraceInfo.rawIntsAtPeak = RayTraceInfo.rawIntsAtPeak(indFiltered);
                    RayTraceInfo.fitIntsAtPeak = RayTraceInfo.fitIntsAtPeak(indFiltered);
                    RayTraceInfo.fitOffset = RayTraceInfo.fitOffset(indFiltered);
                else
                    InterfaceWidthFromFit = InterfaceWidthFit;
                end
            end
        end
        function [recon] = reconCoords(recon)
            [coordsRec,recon.errCoordsRec,recon.fitWidth,~,recon.rayTraceInfo] = ...
                DropRecon3Dv1.reconstructCoords(recon.rayStart,recon.traceVectors,recon.fInterp,...
                recon.fitAlongNormals,recon.roughHEst,recon.nPointsPerPixel,...
                recon.IQRcoordFilter,recon.medianFilterSizeInPixels);
            
            recon.x = coordsRec(:,1);
            recon.y = coordsRec(:,2);
            recon.z = coordsRec(:,3);
            if length(recon.z)>9
                [recon.center,~,~,~,~] = DropRecon3Dv1.ellipsoid_fit( coordsRec );
                recon.xCntrd = recon.x-recon.center(1);
                recon.yCntrd = recon.y-recon.center(2);
                recon.zCntrd = recon.z-recon.center(3);
                recon.message = 'Successful reconstruction of coordinates!';
            else
                recon.message = ['Fewer than 9 coordinates were ',...
                    'reconstructed. Unable to determine elliptical fit.'];
                return
            end
        end
        function [fittedCoords, errFittedCoords,rayTraceInfo,F,sd_est,...
                H_out] = getFittedCoordsCPP(RStart,RV,NPointsPerPixel,F,...
                ts1,ts2,rVscale1,rVscale2,H_est)
            
            %
            NTraces=length(RV);
            fitIntensities = cell(NTraces,1);
            rawIntensities = cell(NTraces,1);
            residuals = cell(NTraces,1);
            fittedCoords = zeros(NTraces,3);
            fWidths = zeros(NTraces,1);
            H_out = zeros(NTraces,1);
            errFittedWidth = zeros(NTraces,1);
            errFittedCenter = zeros(NTraces,1);
            errFittedAmp = zeros(NTraces,1);
            errFittedOffset = zeros(NTraces,1);
            fitIntsAtPeak = zeros(NTraces,1);
            rawIntsAtPeak = zeros(NTraces,1);
            fitOffset = zeros(NTraces,1);
            if NTraces>2.3e4 && license('test', 'Distrib_Computing_Toolbox')
                
                parfor i=1:NTraces
                    rStart1=RStart(i,:);
                    rV0=RV(i,:);
                    traceResult = DropRecon3Dv1.traceRaysCPP(i,rStart1,rV0,...
                        NPointsPerPixel,F,ts1,ts2,rVscale1,rVscale2,H_est);
                    fittedCoords(i,:) = traceResult.fCoord;
                    fWidths(i) = traceResult.fWidth;
                    H_out(i) = traceResult.H_out;
                    errFittedWidth(i) = traceResult.errFittedWidth;
                    errFittedCenter(i) = traceResult.errFittedCenter;
                    errFittedAmp(i) = traceResult.errFittedAmp;
                    errFittedOffset(i) = traceResult.errFittedOffset;
                    fitIntensities{i} = traceResult.fitIntensities;
                    rawIntensities{i} = traceResult.rawIntensities;
                    residuals{i} = traceResult.residuals;
                    fitIntsAtPeak(i) = traceResult.fitIntAtPeak;
                    rawIntsAtPeak(i) = traceResult.rawIntAtPeak;
                    fitOffset(i) = traceResult.fitOffset;
                    
                end
            else
                for i=1:NTraces
                    rStart1=RStart(i,:);
                    rV0=RV(i,:);
                    traceResult = DropRecon3Dv1.traceRaysCPP(i,rStart1,rV0,...
                        NPointsPerPixel,F,ts1,ts2,rVscale1,rVscale2,H_est);
                    fittedCoords(i,:) = traceResult.fCoord;
                    fWidths(i) = traceResult.fWidth;
                    H_out(i) = traceResult.H_out;
                    errFittedWidth(i) = traceResult.errFittedWidth;
                    errFittedCenter(i) = traceResult.errFittedCenter;
                    errFittedAmp(i) = traceResult.errFittedAmp;
                    errFittedOffset(i) = traceResult.errFittedOffset;
                    fitIntensities{i} = traceResult.fitIntensities;
                    rawIntensities{i} = traceResult.rawIntensities;
                    residuals{i} = traceResult.residuals;
                    fitIntsAtPeak(i) = traceResult.fitIntAtPeak;
                    rawIntsAtPeak(i) = traceResult.rawIntAtPeak;
                    fitOffset(i) = traceResult.fitOffset;
                end
                
            end
            fitIntensities = [fitIntensities{:}];
            rawIntensities = [rawIntensities{:}];
            residuals = [residuals{:}];
            
            rayTraceInfo.fitIntensities = fitIntensities(:);
            rayTraceInfo.rawIntensities = rawIntensities(:);
            rayTraceInfo.residuals = residuals(:);
            
            sd_est = fWidths;
            
            errFittedCoords.width=errFittedWidth;
            errFittedCoords.center=errFittedCenter;
            errFittedCoords.amp=errFittedAmp;
            errFittedCoords.offset=errFittedOffset;
            
            rayTraceInfo.fitIntsAtPeak = fitIntsAtPeak;
            rayTraceInfo.rawIntsAtPeak = rawIntsAtPeak;
            rayTraceInfo.fitOffset = fitOffset;
            %disp('Time for tracing rays:')
        end
        function [traceResult] = traceRaysCPP(i,rStart1,rV0,...
                NPointsPerPixel,F,ts1,ts2,rVscale1,rVscale2,H_est)
            
            traceResult.fCoord = [nan, nan, nan];
            traceResult.fWidth = nan;
            traceResult.H_out = nan;
            traceResult.errFittedWidth = nan;
            traceResult.errFittedCenter = nan;
            traceResult.errFittedAmp = nan;
            traceResult.errFittedOffset = nan;
            traceResult.fitIntensities = nan;
            traceResult.rawIntensities = nan;
            traceResult.residuals = nan;
            traceResult.fitIntAtPeak = nan;
            traceResult.rawIntAtPeak = nan;
            traceResult.fitOffset = nan;
            
            %trace line radially from center
            
            %are start and direction vectors useful?
            if (sum(isnan(rStart1))==0 && sum(isnan(rV0))==0)
                %%
                %step 1: prepares the data to be fitted, find coarse guess of
                %center and width
                
                [X,Y,rStart2] = DropRecon3Dv1.prepFitData(rStart1,rV0,...
                    NPointsPerPixel,F,ts1,ts2,rVscale1,rVscale2);
                if ~isempty(nonzeros(X<0))
                    warning('X has no elements. Empty Ray.')
                end
                
                %startP: sigma mu amp offset of Gaussian
                count = 0;
                while length(X)<=10
                    % X is too small, increasing ray length and linear sampling density
                    if count == 0;
                        NPPP = NPointsPerPixel;
                        rVs2 = rVscale2;
                    end
                    NPPP = NPPP +1;
                    rVs2 = rVs2+1;
                    [X,Y,rStart2] = DropRecon3Dv1.prepFitData(rStart1,rV0,...
                        NPPP,F,ts1,ts2,rVscale1,rVs2);
                    count = count + 1;
                    if count == 10
                        break
                    end
                end
                [startP] = DropRecon3Dv1.getInitialFpar(X,Y);
                boolKeep = ~isnan(Y);
                Y=Y(boolKeep);X=X(boolKeep);
                usefulPoints=length(find(isnan(Y)==0  ));
                if usefulPoints>10 && sum(isnan(startP))==0
                    %%
                    %step 2: fits 1D Gaussians using gsl library for speed
                    [fOut,errfOut, rchisquare]=fitGauss(startP,X,Y);
                    fWidth=fOut(1);
                    fCenter=fOut(2);
                    fAmp=fOut(3);
                    fOffset=fOut(4);
                    errfWidth = errfOut(1);
                    errfCenter=errfOut(2);
                    errfAmp = errfOut(3);
                    errfOffset = errfOut(4);
                    YFit=fAmp*exp(-1/2*(X-fCenter).^2/fWidth.^2)+fOffset;
                    %calculate absolute coordinate
                    unit_rV0 = rV0./sqrt(sum(rV0.^2));
                    fitCoordinate = rStart2 + fCenter*unit_rV0;
                    if fitCoordinate(1)>min(F.GridVectors{1})&&...
                            fitCoordinate(1)<max(F.GridVectors{1})&&...
                            fitCoordinate(2)>min(F.GridVectors{2})&&...
                            fitCoordinate(2)<max(F.GridVectors{2})&&...
                            fitCoordinate(3)>min(F.GridVectors{3})&&...
                            fitCoordinate(3)<max(F.GridVectors{3})
                        traceResult.fCoord = fitCoordinate;
                        traceResult.fWidth = fWidth;
                        traceResult.H_out = H_est(i);
                        traceResult.errFittedWidth = errfWidth;
                        traceResult.errFittedCenter = errfCenter;
                        traceResult.errFittedAmp = errfAmp;
                        traceResult.errFittedOffset = errfOffset;
                        traceResult.fitIntensities = YFit;
                        traceResult.rawIntensities = Y;
                        traceResult.residuals = (YFit-Y);
                        traceResult.fitIntAtPeak = fAmp+fOffset;
                        indPeak = find(YFit==max(YFit),1);
                        traceResult.rawIntAtPeak = Y(indPeak);
                        traceResult.fitOffset = fOffset;
                    end
                end
            end
        end
        function [X2,Y2,rStart2] = prepFitData(rStart1,rV0,...
                NPointsPerPixel,F,ts1,ts2,rVscale1,rVscale2)
            %% step 1: trace initial line to find peak center and width
            %get trace vector
            rV1=rVscale1*rV0;
            %generate Coords of Line
            [Pq1, lv1, ldv1,NPoints1,indKeep1]...
                =DropRecon3Dv1.generateLineCoords(rStart1,rV1,NPointsPerPixel,F);
            X1=linspace(0,lv1,NPoints1);
            if length(indKeep1)~=length(X1)
                disp(indKeep1)
                disp(X1)
                %delete(gcp)
                error('length(indPos1) is different from length(X1)')
            end
            X1=X1(indKeep1);
            %interpolates intensities along line
            Y1=F(Pq1(:,indKeep1)');
            %eliminate negative values
            indNeg1= Y1<0;
            %indMaxAmp=find(Y1>255);
            Y1(indNeg1)=nan;
            amp0=max(Y1);
            %find first guess for peak center
            indCenter0=round(mean(find(Y1>ts1*amp0)));
            %if ~isemty(indCenter0)
            if ~isnan(indCenter0)&& indCenter0>0
                center0=X1( indCenter0);
            else
                center0=0;
            end
            %find first guess for peak width
            indCenter0std=round(std(find(Y1>ts1*amp0)));
            center0std=ldv1*indCenter0std;
            %refine peak center by looking at its surrounding
            %take only values within 1 stdev of peak
            indForCenter1=find(abs(X1-center0)<1*center0std);
            amp0b=max(Y1(indForCenter1));
            indCenter1=round(mean(find(Y1(indForCenter1)>ts2*amp0b)));
            if ~isnan(indCenter1)
                center1=X1( indForCenter1(indCenter1));
            else
                center1=0;
            end
            %% step 2: trace line which is centered around the peak determined in step 1
            %absolute coordinate of the peak
            interface=rStart1+center1*1/lv1*rV1;
            %trace line with interface centered
            rV2=rVscale2*rV0;
            rStart2=interface-1/2*rV2;
            %             rEnd2=rStart2+rV2;
            %             while rEnd2(1)<min(F.GridVectors{1})
            %                 %pause(1)
            %                 rV2=rV2*0.9;
            %                 rStart2=interface-1/2*rV2;
            %                 rEnd2=rStart2+rV2;
            %             end
            %             while rEnd2(2)<min(F.GridVectors{2})
            %                % pause()
            %                rV2=rV2*0.9;
            %                rStart2=interface-1/2*rV2;
            %                rEnd2=rStart2+rV2;
            %             end
            %             while rEnd2(3)<min(F.GridVectors{3})
            %                % pause()
            %                rV2=rV2*0.9;
            %                rStart2=interface-1/2*rV2;
            %                rEnd2=rStart2+rV2;
            %             end
            
            [Pq2, lv2, ldv2,NPoints2,indKeep2]=...
                DropRecon3Dv1.generateLineCoords(rStart2,rV2,NPointsPerPixel,F);
            X2 = linspace(0,lv2,length(Pq2));
            X2=X2(indKeep2);
            if length(X2)<3
                warning('ray is too small')
            end
            Y2=F(Pq2');
            Y2=Y2(indKeep2);
            Y2=Y2';
            indNegY2= Y2<0;
            indNegX2= X2<0;
            Y2(indNegY2)=nan;
            Y2(indNegX2)=nan;
        end
        function [lCoords,lv,ldv,NPoints,indKeep] = ...
                generateLineCoords( rStart, rV, NPointsPerPixel,F)
            
            rVLengthEst=norm(rV);
            
            drV = rV./(norm(rV)*NPointsPerPixel);
            
            drVLength=norm(drV);
            NPoints = ceil(rVLengthEst/drVLength);
            rVLength = drVLength*NPoints;
            rq=ones(3,NPoints);
            indNeg=[];
            
            for i=1:3
                %select directions where vector does not change
                if abs(drV(i))<10^(-6)
                    rq(i,:)=(rStart(i)*ones(NPoints,1));
                else
                    rq(i,:)=linspace(0,drV(i)*NPoints,NPoints)+rStart(i);
                end
            end
            
            %checks x y z to be positive
            indPos=sum(rq<0)==0;
            %checks x y z to be within limits
            CoordsMax = [max(F.GridVectors{1});max(F.GridVectors{2});max(F.GridVectors{3})];
            indAcceptable = (sum([rq(1,:)>CoordsMax(1);rq(2,:)>CoordsMax(2);rq(3,:)>CoordsMax(3)],1)==0).*(indPos);
            indKeep = logical(indAcceptable);
            
            lCoords=rq;
            lv=rVLength;
            ldv=drVLength;
            %indPos=indP;
        end
        function [startP] = getInitialFpar(X2,Y2)
            %startP: sigma mu amp offset of Gaussian
            if length(X2) < 3
                warning('X2 is too short; the ray is way too short.')
                startP = NaN;
                return
            end
            amp2=max(Y2);
            dX2=(X2(2)-X2(1));
            %calc mean, stdev of whole trace first guess of peak center and width
            indCenter2=round(mean(find(Y2>0.4*amp2)));
            indCenter2std=round(std(find(Y2>0.1*amp2)));
            
            if isempty(indCenter2) ||isnan(indCenter2)
                startP=[nan nan nan nan];
            else
                center2=X2( indCenter2);
                center2std=dX2*indCenter2std;
                %find coords with intens. values around the peak
                indForCenter2=find(abs(X2-center2)<1.5*center2std);
                indCenter2=round(mean(find(Y2(indForCenter2)>0.4*amp2)));
                
                if isnan(indCenter2)
                    startP=[nan nan nan nan];
                else
                    center2=X2( indForCenter2(indCenter2));
                    width=std((find(Y2>0.2*amp2)))*dX2;
                    
                    %if width > 2*sdmax
                    %width = sdmax;
                    %    warning(['Width of interface seems to be more than',...
                    %        ' twice expected value.'])
                    %end
                    offset=0;
                    %plot(gfitX2,gfitYRaw2)
                    %order of parameters to be fed into CPP startP=[sigma mu amp offset]
                    startP=[width center2 amp2 offset];
                end
            end
        end
        function [indFiltered] = filterInterquartileRange(RawVal,accptRange)
            median(RawVal);
            RawValNotNan=RawVal(~isnan(RawVal)&(RawVal~=Inf));
            
            % compute 25th percentile (first quartile)
            Q(1) = median(RawVal(RawVal<median(RawValNotNan)));
            
            % compute 50th percentile (second quartile)
            Q(2) = median(RawValNotNan);
            
            % compute 75th percentile (third quartile)
            Q(3) = median(RawVal(RawVal>median(RawValNotNan)));
            
            % compute Interquartile Range (IQR)
            IQR = Q(3)-Q(1);
            
            %dataInRange=((Q(2)-accptRange*Q(3))<RawVal) & (RawVal < (Q(2)+accptRange*Q(3)));
            dataInRange=((Q(1)-accptRange*IQR)<RawVal) & (RawVal < (Q(3)+accptRange*IQR));
            length(find(dataInRange));
            length(RawVal);
            
            indFiltered=dataInRange;
        end
        function [NormalsOut,EigenvectorsOut,EigenValuesOut] = ...
                getNormals(CoordsInAll,CoordsInSubset,pKDtree,kNn,...
                patchRadius)
            
            NS_kd = pKDtree;
            %indnnPoints=[];
            %coms=zeros(length(CoordsInSubset),3);
            %patchRadius=[];
            Eigenvalues=zeros(length(CoordsInSubset),3);
            Eigenvectors=zeros(length(CoordsInSubset),3,3);
            Normals=zeros(length(CoordsInSubset),3);
            
            while length(CoordsInSubset) < kNn
                kNn = kNn - 1;
                warning(['Too few points to find normals. Reducing ',...
                    'number of nearest neighbor search by 1.'])
            end
            
            if nargin < 5
                [IDX,DISTS] = knnsearch(NS_kd,CoordsInSubset,'k',kNn);
            else
                [IDX,DISTS] = rangesearch(NS_kd,CoordsInSubset,patchRadius);
            end
            
            CoordsInNear = cell(length(CoordsInSubset),1);
            
            for i = 1:length(CoordsInSubset)
                idx = IDX(i,:);
                dists = DISTS(i,:);
                if nargin > 4
                    dists = dists{1};
                    idx = idx{1};
                end
                indnnPoints=idx;
                indnnPoints=indnnPoints(dists~=0);
                CoordsInNear{i} = CoordsInAll(indnnPoints,:);
            end
            %
            ArgNumber = nargin;
            numCoords = length(CoordsInSubset);
            if numCoords>1.5e5 && license('test', 'Distrib_Computing_Toolbox')
                parfor i=1:length(CoordsInSubset)
                    %search nearest neighbors
                    idx = IDX(i,:);
                    dists = DISTS(i,:);
                    if ArgNumber > 7
                        dists = dists{1};
                        idx = idx{1};
                    end
                    indnnPoints=idx;
                    indnnPoints=indnnPoints(dists~=0);
                    nnPoints=CoordsInNear{i};
                    nNDist=nonzeros(dists);
                    usefulPoints=length(nnPoints(:,1));
                    if (usefulPoints<6)
                        Normals(i,:)=[nan nan nan];
                    else
                        [ eigenvectors,eigenvalues] ...
                            = DropRecon3Dv1.getPrincipalAxes(CoordsInSubset(i,:),...
                            nnPoints);
                        Eigenvectors(i,:,:)=eigenvectors;
                        Eigenvalues(i,:)=eigenvalues;
                        Normals(i,:)=eigenvectors(:,1)';
                    end
                end
                %delete(poolobj)
                %
            else
                %
                for i=1:length(CoordsInSubset)
                    %search nearest neighbors
                    idx = IDX(i,:);
                    dists = DISTS(i,:);
                    if ArgNumber > 7
                        dists = dists{1};
                        idx = idx{1};
                    end
                    indnnPoints=idx;
                    indnnPoints=indnnPoints(dists~=0);
                    nnPoints=CoordsInNear{i};
                    nNDist=nonzeros(dists);
                    usefulPoints=length(nnPoints(:,1));
                    if (usefulPoints<6)
                        Normals(i,:)=[nan nan nan];
                    else
                        [ eigenvectors,eigenvalues] ...
                            = DropRecon3Dv1.getPrincipalAxes(CoordsInSubset(i,:),...
                            nnPoints);
                        Eigenvectors(i,:,:)=eigenvectors;
                        Eigenvalues(i,:)=eigenvalues;
                        Normals(i,:)=eigenvectors(:,1)';
                    end
                end
                %
            end
            
            EigenvectorsOut=Eigenvectors;
            EigenValuesOut=Eigenvalues;
            NormalsOut=Normals;
        end
        function [eigenvectors,eigenvalues] = getPrincipalAxes(qPoint,nnPoints)
            %%  obtain principal axes
            
            usefulPoints=length(nnPoints(:,1));
            COM=mean(nnPoints);
            
            nofnnPoints=length(nnPoints(:,1));
            nnPointsqPoint=nnPoints-ones(nofnnPoints,1)*qPoint;
            d=sqrt(nnPointsqPoint(:,1).^2 + nnPointsqPoint(:,2).^2 + ...
                nnPointsqPoint(:,3).^2);
            
            h2=max(d)^2;
            weightr=exp((-1./2).*(d.^2)./h2);
            nnPointsWM = (weightr*ones(1,3)).*(nnPoints - ...
                ones(usefulPoints,1)*COM);
            
            [~, eigenvectors,eigenvalues]=DropRecon3Dv1.getEigen( nnPointsWM );
        end
        function [normal,eigenvectors,eigenvalues] = getEigen( PointsIn )
            
            RR=zeros(3,3);
            
            % This FOR loop runs in serial because GETEIGEN is called by a
            % worker in a parfor loop
            for l=1:length(PointsIn(:,1))
                rr=PointsIn(l,:)'*PointsIn(l,:);
                RR=RR+rr;
            end
            %eigenvectors and eigenvalues
            %check for nan and inf in matrix
            if (sum(sum(~isnan(RR)))==9 && sum(sum(~isinf(RR)))==9)
                [V,D] = eig(RR);
                eigenvalues=diag(D);
                lev=eigenvalues(1);
                mev=eigenvalues(2);
                %if mev/lev>1.1
                normal=V(:,1)';
                normal=1/norm(normal).*normal;
                eigenvectors=V;
            else
                normal=[nan,nan,nan];
                eigenvectors=[nan,nan,nan;nan,nan,nan;nan,nan,nan];
                eigenvalues=[nan,nan,nan];
            end
        end
        function [NormalsOut] = setNormalOrientationRV( NormalsIn,RVIn )
            % Last edit: 19 Aug 2016
            numNormals = length(NormalsIn(:,1));
            NOut = nan(numNormals,3);
            if numNormals>1.5e6 && license('test', 'Distrib_Computing_Toolbox')
                parfor k=1:numNormals
                    orientation=sign(dot(NormalsIn(k,:),RVIn(k,:)));
                    NOut(k,:)=orientation*NormalsIn(k,:);
                end
                %delete(poolobj)
                %
            else
                %
                for k=1:numNormals
                    orientation=sign(dot(NormalsIn(k,:),RVIn(k,:)));
                    NOut(k,:)=orientation*NormalsIn(k,:);
                end
                %
            end
            NormalsOut=NOut;
        end
        function [RayStarts,Rays] = getRays(rayLength,RayDirections,CoordsInterface)
            [nRays, nCols] = size(RayDirections);
            Rays = nan(nRays, nCols);
            RayStarts = nan(nRays, nCols);
            if nRays>2e6 && license('test', 'Distrib_Computing_Toolbox')
                parfor i=1:length(RayDirections)
                    Rays(i,:)=RayDirections(i,:);
                    Rays(i,:)=rayLength*1/norm(Rays(i,:))*Rays (i,:);
                    RayStarts(i,:)=CoordsInterface(i,:)-1/2.*Rays (i,:);
                end
                %delete(poolobj)
                %
            else
                %
                for i=1:length(RayDirections)
                    Rays(i,:)=RayDirections(i,:);
                    Rays(i,:)=rayLength*1/norm(Rays(i,:))*Rays (i,:);
                    RayStarts(i,:)=CoordsInterface(i,:)-1/2.*Rays (i,:);
                end
                %
            end
        end
        function [CoordsOut] = rotateCoordinates(CoordsIn,PolarAngle,AzimuthalAngle)
            CoordsX_In = CoordsIn(:,1);
            CoordsY_In = CoordsIn(:,2);
            CoordsZ_In = CoordsIn(:,3);
            
            centerX = mean( CoordsX_In );
            centerY = mean( CoordsY_In );
            centerZ = mean( CoordsZ_In );
            
            COM = [centerX,centerY,centerZ];
            
            CenteredCoords = CoordsIn-ones(length(CoordsX_In(:,1)),1)*COM;
            CenteredCoordsX = CenteredCoords(:,1);
            CenteredCoordsY = CenteredCoords(:,2);
            CenteredCoordsZ = CenteredCoords(:,3);
            
            Rot_z = [cos(AzimuthalAngle),-sin(AzimuthalAngle),                 0;...
                sin(AzimuthalAngle), cos(AzimuthalAngle),                 0;...
                0,                    0,                1];
            
            Rot_y = [    cos(PolarAngle),                    0,  sin(PolarAngle);...
                0,                    1,                 0;...
                -sin(PolarAngle),                    0, cos(PolarAngle)];
            
            R_gen = Rot_y*Rot_z;
            
            CenteredCoordsX_rot = R_gen(1,1)*(CenteredCoordsX)+R_gen(1,2)*(CenteredCoordsY)+R_gen(1,3)*(CenteredCoordsZ);
            CenteredCoordsY_rot = R_gen(2,1)*(CenteredCoordsX)+R_gen(2,2)*(CenteredCoordsY)+R_gen(2,3)*(CenteredCoordsZ);
            CenteredCoordsZ_rot = R_gen(3,1)*(CenteredCoordsX)+R_gen(3,2)*(CenteredCoordsY)+R_gen(3,3)*(CenteredCoordsZ);
            
            CoordsX_Out = CenteredCoordsX_rot+centerX;
            CoordsY_Out = CenteredCoordsY_rot+centerY;
            CoordsZ_Out = CenteredCoordsZ_rot+centerZ;
            CoordsOut = [CoordsX_Out,CoordsY_Out,CoordsZ_Out];
        end
        function [center,radii,evecs,v,chi2] = ellipsoid_fit(coords)
            % Fits ellipsoid to a set of coordinates
            %
            %   Adapted from ellipsoid_fit_new, written by Yury Petrov,
            %       Oculus VR, September, 2015
            %
            % Output:
            % * center    -  ellispoid center coordinates [xc; yc; zc]
            % * radii     -  ellipsoid radii [a; b; c]
            % * evecs     -  ellipsoid radii directions as columns of the 3x3 matrix
            % * v         -  the 10 parameters describing the ellipsoid algebraically:
            %                Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
            % * chi2      -  residual sum of squared errors (chi^2), this chi2 is in the
            %                coordinate frame in which the ellipsoid is a unit sphere.
            %
            % Adapted by: Elijah Shelton
            % Last edit 21 January, 2016
            %
            narginchk( 1, 3 ) ;  % check input arguments
            if nargin == 1
                equals = ''; % no constraints by default
            end
            if size(coords,2) ~= 3
                error( 'COORDS must have three columns' );
            else
                x = coords(:,1);
                y = coords(:,2);
                z = coords(:,3);
            end
            
            % check that COORDS contains at least 9 data points
            if length( x ) < 9 && strcmp( equals, '' )
                error( 'At least 9 points are required to fit a unique ellipsoid' );
            end
            
            % fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
            % 2Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
            % parameter
            if strcmp( equals, '' )
                D = [ x .* x + y .* y - 2 * z .* z, ...
                    x .* x + z .* z - 2 * y .* y, ...
                    2 * x .* y, ...
                    2 * x .* z, ...
                    2 * y .* z, ...
                    2 * x, ...
                    2 * y, ...
                    2 * z, ...
                    1 + 0 * x ];  % ndatapoints x 9 ellipsoid parameters
            end
            % solve the normal system of equations
            d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
            u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations
            if strcmp( equals, '' )
                v(1) = u(1) +     u(2) - 1;
                v(2) = u(1) - 2 * u(2) - 1;
                v(3) = u(2) - 2 * u(1) - 1;
                v( 4 : 10 ) = u( 3 : 9 );
            end
            v = v';
            
            % form the algebraic form of the ellipsoid
            A = [ v(1) v(4) v(5) v(7); ...
                v(4) v(2) v(6) v(8); ...
                v(5) v(6) v(3) v(9); ...
                v(7) v(8) v(9) v(10) ];
            % find the center of the ellipsoid
            center = -A( 1:3, 1:3 ) \ v( 7:9 );
            % form the corresponding translation matrix
            T = eye( 4 );
            T( 4, 1:3 ) = center';
            % translate to the center
            R = T * A * T';
            % solve the eigenproblem
            [ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
            radii = sqrt( 1 ./ diag( evals ) );
            
            % calculate difference of the fitted points from the actual data normalized by the ellipsoid radii
            d = [ x - center(1), y - center(2), z - center(3) ]; % shift data to origin
            d = d * evecs; % rotate to cardinal axes of the ellipsoid;
            d = [ d(:,1) / radii(1), d(:,2) / radii(2), d(:,3) / radii(3) ]; % normalize to the ellipsoid radii
            chi2 = sum( abs( 1 - sum( d.^2, 2 ) ) );
            
            if abs( v(end) ) > 1e-6
                v = -v / v(end); % normalize to the more conventional form with constant term = -1
            else
                v = -sign( v(end) ) * v;
            end
            
        end
        function [recon] = findAxiallySymmetricEllipsoid(recon)
            % Inspired by
            % ellipsoid_fit_new // by Yury Petrov, Oculus VR // Date: September, 2015
            % FIRST, fit a tri-axial ellipsoid
            x = recon.x;
            y = recon.y;
            z = recon.z;
            D = [ x .* x + y .* y - 2 * z .* z, ...
                x .* x + z .* z - 2 * y .* y, ...
                2 * x .* y, ...
                2 * x .* z, ...
                2 * y .* z, ...
                2 * x, ...
                2 * y, ...
                2 * z, ...
                1 + 0 * x ];  % ndatapoints x 9 ellipsoid parameters
            % solve the normal system of equations
            d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
            u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations
            % find the ellipsoid parameters
            % convert back to the conventional algebraic form
            v = zeros(10,1);
            v(1) = u(1) +     u(2) - 1;
            v(2) = u(1) - 2 * u(2) - 1;
            v(3) = u(2) - 2 * u(1) - 1;
            v( 4 : 10 ) = u( 3 : 9 );
            % form the algebraic form of the ellipsoid
            A = [ v(1) v(4) v(5) v(7); ...
                v(4) v(2) v(6) v(8); ...
                v(5) v(6) v(3) v(9); ...
                v(7) v(8) v(9) v(10) ];
            % find the center of the ellipsoid
            recon.center = -A( 1:3, 1:3 ) \ v( 7:9 );
            % form the corresponding translation matrix
            T = eye( 4 );
            T( 4, 1:3 ) = recon.center';
            % translate to the center
            R = T * A * T';
            % solve the eigenproblem
            [ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
            radii = sqrt( 1 ./ diag( evals ) );
            %% SECOND, determine which radii is most different from the
            % other two, and rotate the coordinates into this frame.
            B = (radii-mean(radii)).^2;
            indPA = find(B==max(B));
            [aziAxialEllipsoid, eleAxialEllipsoid, ~] = cart2sph(evecs(1,indPA),...
                evecs(2,indPA),evecs(3,indPA));
            polAxialEllipsoid = pi/2-eleAxialEllipsoid;
            
            recon.aziAxialEllipsoid = aziAxialEllipsoid;
            recon.eleAxialEllipsoid = eleAxialEllipsoid;
            recon.polAxialEllipsoid = polAxialEllipsoid;
            
            recon.xCntrd = recon.x-recon.center(1);
            recon.yCntrd = recon.y-recon.center(2);
            recon.zCntrd = recon.z-recon.center(3);
            coordsLocalFrame = DropRecon3Dv1.rotateCoordinates([recon.xCntrd,recon.yCntrd,recon.zCntrd],-polAxialEllipsoid,-aziAxialEllipsoid);
            recon.xLocal = coordsLocalFrame(:,1);
            recon.yLocal = coordsLocalFrame(:,2);
            recon.zLocal = coordsLocalFrame(:,3);
            
            x = recon.xLocal;
            y = recon.yLocal;
            z = recon.zLocal;
            %% THIRD, fit an axially symmetric ellipsoid to the coordinates.
            D = [ x .* x + y .* y - 2 * z .* z, ...
                2 * x, ...
                2 * y, ...
                2 * z, ...
                1 + 0 * x ];  % ndatapoints x 5 ellipsoid parameters
            
            % solve the normal system of equations
            d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
            u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations
            v = zeros(10,1);
            v(1) = u(1) - 1;
            v(2) = u(1) - 1;
            v(3) = -2 * u(1) - 1;
            v = [ v( 1 : 6 ); u( 2 : 5 ) ];
            % form the algebraic form of the ellipsoid
            A = [ v(1) v(4) v(5) v(7); ...
                v(4) v(2) v(6) v(8); ...
                v(5) v(6) v(3) v(9); ...
                v(7) v(8) v(9) v(10) ];
            % find the center of the ellipsoid
            center = -A( 1:3, 1:3 ) \ v( 7:9 );
            % form the corresponding translation matrix
            T = eye( 4 );
            T( 4, 1:3 ) = center';
            % translate to the center
            R = T * A * T';
            % solve the eigenproblem
            [ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
            radii = sqrt( 1 ./ diag( evals ) );
            %% FOURTH, calculate the curvatures of the fitted ellipsoid
            a = radii((evecs*[1;0;0]).^2==max((evecs*[1;0;0]).^2));
            b = radii((evecs*[0;1;0]).^2==max((evecs*[0;1;0]).^2));
            c = radii((evecs*[0;0;1]).^2==max((evecs*[0;0;1]).^2));
            
            recon.aAxiPixels = a;
            recon.bAxiPixels = b;
            recon.cAxiPixels = c;
            
            recon.aAxiMicrons = recon.aAxiPixels*recon.vsx;
            recon.bAxiMicrons = recon.bAxiPixels*recon.vsx;
            recon.cAxiMicrons = recon.cAxiPixels*recon.vsx;
            
            xo = recon.xLocal;
            yo = recon.yLocal;
            zo = recon.zLocal;
            
            u_ref = atan((a*yo)./(b*xo));
            u = u_ref;
            Q1 = (xo>0)&(yo>0);
            Q2 = (xo<0)&(yo>0);
            Q3 = (xo<0)&(yo<0);
            Q4 = (xo>0)&(yo<0);
            
            u(Q1)=u_ref(Q1);
            u(Q2)=u_ref(Q2)+pi;
            u(Q3)=u_ref(Q3)+pi;
            u(Q4)=u_ref(Q4)+2*pi;
            
            cosu = cos(u);
            t = sqrt(1./((zo./c).^2+(xo./(a*cosu)).^2));
            cosv = zo.*t/c;
            sinv = xo.*t./(a*cosu);
            sinu = yo.*t./(b*sinv);
            
            cos2v = cosv.^2-sinv.^2;
            cos2u = cosu.^2-sinu.^2;
            
            recon.H_AxialEllipsoid = ((1/8)*(a*b*c)*(3*(a^2+b^2)+2*c^2+(a^2+b^2-2*c^2)*...
                cos2v-2*(a^2-b^2)*cos2u.*sinv.^2)./...
                ((a^2*b^2*cosv.^2+c^2*(b^2*cosu.^2+a^2*sinu.^2).*sinv.^2).^(3/2)));
            
            
            recon.K_AxialEllipsoid = real(((a*b*c)^2)./((a*b*cosv).^2+(c^2).*...
                ((b*cosu).^2+(a*sinu).^2).*sinv.^2).^2);
            
            recon.HaAxi = (a*b*c*(b^2 + c^2))/(2*(b*c)^3); % in inverse pixels
            recon.HbAxi = (a*b*c*(a^2 + c^2))/(2*(a*c)^3); % in inverse pixels
            recon.HcAxi = (a*b*c*(a^2 + b^2))/(2*(a*b)^3); % in inverse pixels
            
            recon.deltaHacAxi = recon.HaAxi - recon.HcAxi;
            recon.deltaHabAxi = recon.HaAxi - recon.HbAxi;
            recon.deltaHbcAxi = recon.HbAxi - recon.HcAxi;
            
            if recon.deltaHacAxi>0
                recon.prolateOrOblate = 'oblate';
            elseif recon.deltaHacAxi<0
                recon.prolateOrOblate = 'prolate';
            else
                recon.prolateOrOblate = 'neither';
            end
        end
        function [recon] = findGenEllipsoid(recon)
            %% FIRST, fit ellipsoid to coordinates
            coords = [recon.x,recon.y,recon.z];
            [center, radii, evecs, v, chi2 ] = DropRecon3Dv1.ellipsoid_fit(coords);
            %% SECOND, calulate curvatures at poles and determine
            % orientations of axes
            r1 = radii(1);
            r2 = radii(2);
            r3 = radii(3);
            
            recon.aGenPixels = r1;
            recon.bGenPixels = r2;
            recon.cGenPixels = r3;
            
            recon.aGenMicrons = recon.aGenPixels*recon.vsx;
            recon.bGenMicrons = recon.bGenPixels*recon.vsx;
            recon.cGenMicrons = recon.cGenPixels*recon.vsx;
            
            recon.HaGen = (r1*r2*r3*(r2^2 + r3^2))/(2*(r2*r3)^3);
            recon.HbGen = (r1*r3*r3*(r1^2 + r3^2))/(2*(r1*r3)^3);
            recon.HcGen = (r1*r2*r3*(r1^2 + r2^2))/(2*(r1*r2)^3);
            
            for i = 1:3
                [recon.aziGenEllipsoid(i),recon.eleGenEllipsoid(i),~] = cart2sph(evecs(1,i),evecs(2,i),evecs(3,i));
            end
            
            %% THIRD, rotate coordinates into local frame of the third axis
            % of this ellipsoid (this is the shortest axis). Perform the
            % fit again.
            
            recon.xCntrd = recon.x-center(1);
            recon.yCntrd = recon.y-center(2);
            recon.zCntrd = recon.z-center(3);
            
            polAxis3 = pi/2 - recon.eleGenEllipsoid(3);
            coordsLocalFrame = DropRecon3Dv1.rotateCoordinates([recon.xCntrd,recon.yCntrd,recon.zCntrd],-polAxis3,-recon.aziGenEllipsoid(3));
            recon.xLocal = coordsLocalFrame(:,1);
            recon.yLocal = coordsLocalFrame(:,2);
            recon.zLocal = coordsLocalFrame(:,3);
            
            [center, radii, evecs, v, chi2 ] = DropRecon3Dv1.ellipsoid_fit(coordsLocalFrame);
            
            %% FOURTH, rotate the coordinates into a local frame such that
            % the largest semi-axis is along X, the medium semi-axis is
            % along Y, and the shortest semi-axis is along Z.
            
            [azimuth,~,~] = cart2sph(evecs(1,1),evecs(2,1),evecs(3,1));
            coordsLocalFrame = DropRecon3Dv1.rotateCoordinates([recon.xLocal,recon.yLocal,recon.zLocal],0,-azimuth);
            recon.xLocal = coordsLocalFrame(:,1);
            recon.yLocal = coordsLocalFrame(:,2);
            recon.zLocal = coordsLocalFrame(:,3);
            
            %% FIFTH, compute the analytical curvature on the ellipsoid
            
            a = r1;
            b = r2;
            c = r3;
            
            xo = recon.xLocal;
            yo = recon.yLocal;
            zo = recon.zLocal;
            
            u_ref = atan((a*yo)./(b*xo));
            u = u_ref;
            Q1 = (xo>0)&(yo>0);
            Q2 = (xo<0)&(yo>0);
            Q3 = (xo<0)&(yo<0);
            Q4 = (xo>0)&(yo<0);
            
            u(Q1)=u_ref(Q1);
            u(Q2)=u_ref(Q2)+pi;
            u(Q3)=u_ref(Q3)+pi;
            u(Q4)=u_ref(Q4)+2*pi;
            
            cosu = cos(u);
            t = sqrt(1./((zo./c).^2+(xo./(a*cosu)).^2));
            cosv = zo.*t/c;
            sinv = xo.*t./(a*cosu);
            sinu = yo.*t./(b*sinv);
            
            cos2v = cosv.^2-sinv.^2;
            cos2u = cosu.^2-sinu.^2;
            
            recon.H_GenEllipsoid = (a*b*c*(3*(a^2+b^2)+2*c^2+(a^2+b^2-2*c^2)*...
                cos2v-2*(a^2-b^2)*cos2u.*sinv.^2)./...
                (8*(a^2*b^2*cosv.^2+c^2*(b^2*cosu.^2+a^2*sinu.^2).*sinv.^2).^(3/2)));
            recon.K_GenEllipsoid = real(((a*b*c)^2)./((a*b*cosv).^2+(c^2).*...
                ((b*cosu).^2+(a*sinu).^2).*sinv.^2).^2);
        end
        function [recon] = applyMedianFilter2Recon(recon)
            
            if isempty(recon.medFiltPatch.idx)
                [recon.azi,recon.elv,recon.r] = ...
                    cart2sph(recon.xCntrd,recon.yCntrd,recon.zCntrd);
                CoordsIn = [recon.xCntrd,recon.yCntrd,recon.zCntrd];
                NS_kd = KDTreeSearcher(CoordsIn);
                startIndex = 1;
                endIndex = length(recon.x);
                patchSize = recon.medianFilterSizeInPixels;
                %poolObj = parpool;
                r_out = nan(endIndex,1);
                recon_medFiltPatch_idx = cell(endIndex,1);
                recon_r = recon.r;
                if license('test', 'Distrib_Computing_Toolbox')
                    parfor i = startIndex:endIndex
                        %% GET CIRCULAR PATCH COORDINATES
                        queryPoint = CoordsIn(i,:);
                        % find indices of points in the neighborhood
                        [idx,~]=rangesearch(NS_kd,queryPoint,patchSize);
                        idx = idx{1};
                        idx=nonzeros(idx);
                        recon_medFiltPatch_idx{i} = idx;
                        % this is the neighborhood
                        r_patch = recon_r(idx);
                        r_out(i) = median(r_patch);
                    end
                else
                    for i = startIndex:endIndex
                        %% GET CIRCULAR PATCH COORDINATES
                        queryPoint = CoordsIn(i,:);
                        % find indices of points in the neighborhood
                        [idx,~]=rangesearch(NS_kd,queryPoint,patchSize);
                        idx = idx{1};
                        idx=nonzeros(idx);
                        recon_medFiltPatch_idx{i} = idx;
                        % this is the neighborhood
                        r_patch = recon_r(idx);
                        r_out(i) = median(r_patch);
                    end
                end
                recon.medFiltPatch.idx = recon_medFiltPatch_idx;
                
            else
                startIndex = 1;
                endIndex = length(recon.x);
                r_out = nan(endIndex,1);
                for i = startIndex:endIndex
                    %% GET CIRCULAR PATCH COORDINATES
                    idx = recon.medFiltPatch.idx{i};
                    %this is the neighborhood
                    r_patch = recon.r(idx);
                    r_out(i) = median(r_patch);
                end
            end
            recon.r = r_out;
            [recon.xCntrd,recon.yCntrd,recon.zCntrd] = ...
                sph2cart(recon.azi,recon.elv,recon.r);
            recon.x = recon.xCntrd+recon.center(1);
            recon.y = recon.yCntrd+recon.center(2);
            recon.z = recon.zCntrd+recon.center(3);
            
        end
        function [recon] = makeGriddedSurface(recon)
            %% 22 Sept 2016
            %% Elijah Shelton
            %% Interpolate Surface Coordinates for Surface Plots
            [aziRec,elvRec,rRec] = cart2sph(recon.xCntrd,recon.yCntrd,recon.zCntrd);
            
            method = 'natural';
            n=recon.oversamplingFactor*ceil(pi*max(rRec)/recon.d_lambda);
            
            aziVec = linspace(0,2*pi,2*n+1);
            
            polVec = acos(linspace(-1,1,n));
            [recon.polGrid,recon.aziGrid] = meshgrid(polVec,aziVec);
            recon.elGrid = pi/2 - recon.polGrid;
            recon.aziRecForInterp = [(aziRec-2*pi);aziRec;(aziRec+2*pi)];
            recon.elvRecForInterp  = [elvRec;elvRec;elvRec];
            recon.rRecForInterp  = [rRec;rRec;rRec];
            F_r = scatteredInterpolant(...
                recon.aziRecForInterp,recon.elvRecForInterp,...
                recon.rRecForInterp,method);
            recon.rGrid = F_r(recon.aziGrid,recon.elGrid);
            [recon.xGrid,recon.yGrid,recon.zGrid] = ...
                sph2cart(recon.aziGrid,recon.elGrid,recon.rGrid);
        end
        function [anal] = getAnalytical(recon)
            %% 21 Sept 2016
            %% by Elijah Shelton
            %% Get Analytical Values for Reconstructed Coordinates
            coords = [recon.xCntrd,recon.yCntrd,recon.zCntrd];
            coeffs = recon.coeffs;
            [shape] = sphrHarm.getCurvsFromCoords(coords,coeffs);
            anal.t = shape.t(:);
            anal.r = shape.r(:);
            anal.rt = shape.rt(:);
            anal.rtt = shape.rtt(:);
            anal.ks = shape.ks(:);
            anal.kp = shape.kp(:);
            anal.s = shape.s(:);
            anal.H = shape.H(:);
            
            %% Get Analytical Values for Gridded Coordinates
            xGrid = recon.xGrid;
            yGrid = recon.yGrid;
            zGrid = recon.zGrid;
            coordsGrid = [xGrid(:),yGrid(:),zGrid(:)];
            [rows,cols] = size(xGrid);
            [shape] = sphrHarm.getCurvsFromCoords(coordsGrid,coeffs);
            anal.tGrid = reshape(shape.t, [rows cols]);
            anal.rGrid = reshape(shape.r, [rows cols]);
            anal.rtGrid = reshape(shape.rt, [rows cols]);
            anal.rttGrid = reshape(shape.rtt, [rows cols]);
            anal.ksGrid = reshape(shape.ks, [rows cols]);
            anal.kpGrid = reshape(shape.kp, [rows cols]);
            anal.sGrid = reshape(shape.s, [rows cols]);
            anal.HGrid = reshape(shape.H, [rows cols]);
            
            
        end
        function [evalMetric,recon] = evalRecon(recon,anal)
            %% 21 Sept 2016
            %% by Elijah Shelton
            
            evalMetric.rRelErr = (recon.r - anal.r)./recon.r;
            evalMetric.rErr = (recon.r - anal.r);
            evalMetric.sigRelErr = (recon.fitWidth - recon.sig)./recon.sig;
            evalMetric.sigErr = (recon.fitWidth - recon.sig);
            
            evalMetric.rRelErr_mean = mean(evalMetric.rRelErr);
            evalMetric.rErr_mean = mean(evalMetric.rErr);
            evalMetric.sigRelErr_mean = mean(evalMetric.sigRelErr);
            evalMetric.sigErr_mean = mean(evalMetric.sigErr);
            
            evalMetric.rRelErr_std = std(evalMetric.rRelErr);
            evalMetric.rErr_std = std(evalMetric.rErr);
            evalMetric.sigRelErr_std = std(evalMetric.sigRelErr);
            evalMetric.sigErr_std = std(evalMetric.sigErr);
            
            
            recon.r_analGrid = anal.rGrid;
            
        end
        function [measurement] = recon2Measurement(recon)
            fType = 'eps';
            switch recon.message
                case ['Fewer than 9 coordinates were ',...
                        'reconstructed. Unable to determine elliptical fit.']
                    return
                case 'Successful reconstruction of coordinates!';
                    measurement = recon;
                    clear('recon')
                    close all
                    if measurement.findCurvs
                        measurement.numberOfPoints = length(measurement.x);
                        i = 1;
                        measurement.iterNum=i;
                        switch measurement.degPoly
                            case 2
                                measurement=DropRecon3Dv1.circularPatchFitsPoly2(...
                                    measurement);
                            case 3
                                measurement=DropRecon3Dv1.circularPatchFitsPoly3(...
                                    measurement);
                        end
                        measurement = ...
                            DropRecon3Dv1.applyMedianFilter2Measurement(measurement);
                        measurement.patchHistA=measurement.a_patch;
                        measurement.patchHistB=measurement.b_patch;
                        FOM = measurement.GOF.p;
                        measurement.FOM_hist=FOM;
                        noi = measurement.numOfIterations;
                        measurement.continue = 1;
                        if noi<2
                            measurement.continue = 0;
                        end
                        measurement.HRecForInterp = [measurement.H_rec;...
                            measurement.H_rec;measurement.H_rec];
                        F_h = scatteredInterpolant(...
                            measurement.aziRecForInterp,...
                            measurement.elvRecForInterp,...
                            measurement.HRecForInterp,'natural');
                        measurement.HGrid = F_h(...
                            measurement.aziGrid,...
                            measurement.elGrid);
                        suffix = ['Iter',num2str(i)];
                        xGrid = measurement.xGrid;
                        yGrid = measurement.yGrid;
                        zGrid = measurement.zGrid;
                        HGrid = measurement.HGrid;
                        rGrid = measurement.rGrid;
                        
                        while measurement.continue
                            i = i+1;
                            measurement.iterNum=i;
                            suffix = ['Iter',num2str(i)];
                            measurement.FOM = measurement.GOF.p;
                            measurement.x = ...
                                measurement.coordsRecSurface(:,1);
                            measurement.y = ...
                                measurement.coordsRecSurface(:,2);
                            measurement.z = ...
                                measurement.coordsRecSurface(:,3);
                            
                            switch measurement.degPoly
                                case 2
                                    measurement=DropRecon3Dv1.circularPatchFitsPoly2(...
                                        measurement);
                                case 3
                                    measurement=DropRecon3Dv1.circularPatchFitsPoly3(...
                                        measurement);
                            end
                            measurement = ...
                                DropRecon3Dv1.applyMedianFilter2Measurement(...
                                measurement);
                            measurement.HRecForInterp = [measurement.H_rec;...
                                measurement.H_rec;measurement.H_rec];
                            F_h = scatteredInterpolant(...
                                measurement.aziRecForInterp,...
                                measurement.elvRecForInterp,...
                                measurement.HRecForInterp,'natural');
                            measurement.HGrid = F_h(...
                                measurement.aziGrid,...
                                measurement.elGrid);
                            if i>noi;
                                measurement.continue=0;
                            end
                        end
                    end
            end
        end
        function [measurement] = circularPatchFitsPoly2(measurement)
            % Elijah Shelton
            % Last Edited: 31 Aug 2016
            a1 = 0.198737952524216;
            b1 = 0.625861735323673;
            c1 = 1.94278276279535;
            maxHVar = measurement.maxHVar;
            d_lambda = measurement.d_lambda;
            minRp = measurement.minRp;
            maxRp = mean([measurement.aGenPixels,measurement.cGenPixels,...
                measurement.bGenPixels]);
            CoordsIn = [measurement.x,measurement.y,measurement.z];
            if isempty(measurement.H_rec)
                k_est = measurement.H_AxialEllipsoid;
            else
                k1 = measurement.k1;
                k2 = measurement.k2;
                k_est = max([k1';k2'])';
            end
            numberOfPoints = measurement.numberOfPoints;
            
            if length(CoordsIn(:,1))>numberOfPoints
                
                fprintf(['Only up to ', num2str(numberOfPoints),' points ',...
                    'will be considered in the curvature calculation.\n'])
                fprintf(['This dataset contains %d points. Not all',...
                    ' will be considered. The current maximum is set to %d.',...
                    ' Consider selecting ''Custom'' in the dialog box ',...
                    'to increase the number of points to be considered.',...
                    'considered.\n'],length(CoordsIn(:,1)),numberOfPoints)
            end
            numInNest = 50000;
            NS_kd = KDTreeSearcher(CoordsIn);
            center = mean(CoordsIn);
            
            %% GET PATCH ESTIMATES FROM CURVATURE ESTIMATE FROM ECCENTRICITY
            %alpha = std(measurement.rRec-measurement.R_anal);
            if isempty(measurement.H_rec)
                alpha = 1;
            else
                alpha = measurement.Alphas(:,measurement.iterNum-1);
            end
            switch measurement.patchMethod
                case 'SphereRule'
                    patchRadii = c1*(alpha.^a1).*(1./k_est).^b1;
                case 'MinimalRule'
                    patchRadii = minRp*ones(length(k_est),1);
                case 'onePxNormChng'
                    patchRadii = (2./k_est).^(1/2);
            end
            
            patchRadii(patchRadii<minRp)=minRp;
            patchRadii(patchRadii>maxRp)=maxRp;
            PR = patchRadii;
            roughInd = find(isnan(patchRadii));
            numOfRoughPatches = length(roughInd);
            for n = 1:numOfRoughPatches
                i = roughInd(n);
                queryPoint = CoordsIn(i,:);
                [idx,~]=rangesearch(NS_kd,queryPoint,minRp); % this delivers the indices of points in the neighborhood
                idx = idx{1};
                idx=nonzeros(idx);
                PRinPatch = PR(idx);
                PRinPatch = PRinPatch(~isnan(PRinPatch));
                if isempty(PRinPatch)
                    PRinPatch = minRp;
                end
                patchRadii(i) = median(PRinPatch);
            end
            a_patch = nan(numberOfPoints,1);
            b_patch = nan(numberOfPoints,1);
            coordsRecSurface = nan(numberOfPoints,3);
            surfaceNormalsPoly3 = nan(numberOfPoints,3);
            MeanCurv = nan(numberOfPoints,1);
            H_rangeInPatch = nan(numberOfPoints,1);
            MeanCurvStdErr = nan(numberOfPoints,1);
            GausCurv = nan(numberOfPoints,1);
            ChiSqrd = nan(numberOfPoints,1);
            ResidualsInnerMean = nan(numberOfPoints,1);
            ResidualsOuterMean = nan(numberOfPoints,1);
            ResidualsTimesRadius = nan(numberOfPoints,1);
            fittingParameters = nan(numberOfPoints,6);
            LengthScaleVariation = nan(numberOfPoints,1);
            LengthScalePatch = nan(numberOfPoints,1);
            Alphas = nan(numberOfPoints,1);
            Betas = nan(numberOfPoints,1);
            GOF_h = nan(numberOfPoints,1);
            GOF_p = nan(numberOfPoints,1);
            P11F = nan(numberOfPoints,1);
            P20F = nan(numberOfPoints,1);
            P02F = nan(numberOfPoints,1);
            NumPerPatch = nan(numberOfPoints,1);
            numOfNestedLoops = ceil(numberOfPoints/numInNest);
            for n=1:numOfNestedLoops
                startIndex = numInNest*(n-1)+1;
                endIndex = numInNest*n;
                if endIndex>numberOfPoints
                    endIndex = numberOfPoints;
                end
                if license('test', 'Distrib_Computing_Toolbox')
                    parfor i = startIndex:endIndex
                        %% THESE NEED TO BE SLICED
                        queryPoint = CoordsIn(i,:);
                        rp = patchRadii(i);
                        [MeanCurv(i),GausCurv(i),patchRadii(i),...
                            ResidualsInnerMean(i),ResidualsOuterMean(i),...
                            ResidualsTimesRadius(i),H_rangeInPatch(i),...
                            fittingParameters(i,:),a_patch(i),b_patch(i),...
                            coordsRecSurface(i,:),...
                            surfaceNormalsPoly3(i,:),MeanCurvStdErr(i),...
                            ChiSqrd(i),LengthScaleVariation(i),...
                            LengthScalePatch(i),Alphas(i),Betas(i),GOF_h(i),...
                            GOF_p(i),P11F(i),P20F(i),P02F(i),NumPerPatch(i)] = ...
                            DropRecon3Dv1.littleAnalysisLoop(i,CoordsIn,rp,NS_kd,k_est,minRp,maxHVar,center,d_lambda);
                    end
                else
                    for i = startIndex:endIndex
                        
                        %% GET CIRCULAR PATCH COORDINATES
                        queryPoint = CoordsIn(i,:);
                        patchTooBig = 1;
                        rp = patchRadii(i);
                        while patchTooBig
                            [idx,~]=rangesearch(NS_kd,queryPoint,rp); % this delivers the indices of points in the neighborhood
                            patchRadii(i) = rp;
                            idx = idx{1};
                            idx=nonzeros(idx);
                            %this is the neighborhood
                            patchCoordsGF=CoordsIn(idx,:);
                            usefulPoints = length(patchCoordsGF(:,1));
                            H_inPatch = k_est(idx);
                            deltaH_rec = range(H_inPatch);
                            if (rp==minRp)||(usefulPoints==6)
                                patchTooBig = 0;
                                if usefulPoints<6
                                    [idx,dst]=knnsearch(NS_kd,queryPoint,'k',6); % this delivers the indices of points in the neighborhood
                                    idx=nonzeros(idx);
                                    %this is the neighborhood
                                    patchCoordsGF=CoordsIn(idx,:);
                                    patchRadii(i) = max(dst(:));
                                    usefulPoints = length(patchCoordsGF(:,1));
                                end
                            else
                                patchTooBig = (deltaH_rec/k_est(i))>maxHVar;
                                if patchTooBig
                                    rp = rp*(1-2*pi/d_lambda);
                                    if (rp<minRp)||(usefulPoints<6)
                                        rp=minRp;
                                    end
                                end
                            end
                        end
                        if (usefulPoints<6)
                            error('Too few useful points.')
                        end
                        [eigenvectors,~] = DropRecon3Dv1.getPrincipalAxes(...
                            queryPoint,patchCoordsGF);
                        % e1 cooresponds to normal
                        e1=eigenvectors(:,1);e2=eigenvectors(:,2);...
                            e3=eigenvectors(:,3);
                        % fix orientation of normal
                        if sign((center - queryPoint)*e1) == -1
                            e1 = -e1;
                            e2 = -e2;
                            %e3 = e3;
                        end
                        
                        if (sum([isnan(e1)', isnan(e2)', isnan(e3)'])>0)
                            warning('Principal axes are meaningless')
                        end
                        thisCoordinate = [nan, nan, nan];
                        gofH = nan;
                        gofP = nan;
                        p20f = nan;
                        p02f = nan;
                        p11f = nan;
                        N_GF = [nan, nan, nan];
                        usefulPoints = length(patchCoordsGF(:,1));
                        if (usefulPoints<6)
                            warning('Too few points for meaningful fit')
                            H = nan;
                            dH = nan;
                            K = nan;
                            inRes = nan;
                            outRes = nan;
                            RRes = nan;
                            LengthVar = nan;
                            LengthPatch = nan;
                            fPar2 = nan(1,6);
                            chiSq2 = nan;
                            alpha = nan;
                            beta = nan;
                            npp = usefulPoints;
                        else
                            [patchCoordsLF,W_patch,qPoint] ...
                                = DropRecon3Dv1.getCoordsLocalFrame(queryPoint,patchCoordsGF,e1,e2,e3 );
                            %% FIT CIRCULAR PATCH COORDINATES
                            X_patch=patchCoordsLF(:,1);
                            Y_patch=patchCoordsLF(:,2);
                            Z_patch=patchCoordsLF(:,3);
                            
                            [fPar2,errfPar2,chiSq2,residuals2] = ...
                                DropRecon3Dv1.getPolyFit2(X_patch,...
                                Y_patch,Z_patch,W_patch);
                            
                            p20f = fPar2(1);
                            p02f = fPar2(2);
                            p11f = fPar2(3);
                            p10 = fPar2(4);
                            p01 = fPar2(5);
                            p00 = fPar2(6);
                            
                            qx0=qPoint(1);qy0=qPoint(2);
                            % This is how we turn the fit data into mean curvature
                            
                            [H,dH]=DropRecon3Dv1.getMeanCurvaturePoly2(...
                                fPar2,qx0,qy0,errfPar2);
                            % This is how we turn the fit data into gaussian curvature
                            [K]=DropRecon3Dv1.getGaussianCurvaturePoly2(...
                                fPar2,qx0,qy0);
                            
                            [nxLF,nyLF,nzLF]=DropRecon3Dv1.getNormalPoly2(...
                                fPar2,qx0,qy0);
                            N_LF = [nxLF,nyLF,nzLF];
                            
                            alpha = std(residuals2);
                            beta = mean(residuals2);
                            [gofH,gofP] = chi2gof(residuals2);
                            R_patch = sqrt(X_patch.^2+Y_patch.^2+...
                                Z_patch.^2);
                            R_median = median(R_patch);
                            inRes = mean(residuals2(R_patch<R_median));
                            outRes = mean(residuals2(R_patch>R_median));
                            RRes = sum(residuals2.*R_patch)./sum(R_patch);
                            LengthVar = nan;
                            LengthPatch = patchRadii(i);
                            npp = length(X_patch);
                            qzFit = p20f*qx0^2+p02f*qy0^2+...
                                p11f*qx0*qy0+p10*qx0+p01*qy0+p00;
                            qPointLF = [qx0,qy0,qzFit];
                            N_GF = N_LF/[e2,e3,e1];
                            thisCoordinate = qPointLF/[e2,e3,e1]+...
                                mean(patchCoordsGF);
                        end
                        coordsRecSurface(i,:) = thisCoordinate;
                        surfaceNormalsPoly3(i,:) = N_GF;
                        H_rangeInPatch(i) = deltaH_rec;
                        MeanCurv(i) = H;
                        MeanCurvStdErr(i) = dH;
                        GausCurv(i) = K;
                        ChiSqrd(i) = chiSq2;
                        ResidualsInnerMean(i) = inRes;
                        ResidualsOuterMean(i) = outRes;
                        ResidualsTimesRadius(i) = RRes;
                        LengthScaleVariation(i) = LengthVar;
                        LengthScalePatch(i) = LengthPatch;
                        a_patch(i) = LengthPatch;
                        b_patch(i) = LengthPatch;
                        fittingParameters(i,:) = fPar2;
                        Alphas(i) = alpha;
                        Betas(i) = beta;
                        GOF_h(i) = gofH;
                        GOF_p(i) = gofP;
                        P11F(i) = p11f;
                        P02F(i) = p02f;
                        P20F(i) = p20f;
                        NumPerPatch(i) = npp;
                    end
                end
            end
            
            Residuals.innerMean = ResidualsInnerMean;
            Residuals.outerMean = ResidualsOuterMean;
            Residuals.timesR =  ResidualsTimesRadius;
            
            k1 = MeanCurv-sqrt(MeanCurv.^2-GausCurv); % smaller principal curvature
            k2 = MeanCurv+sqrt(MeanCurv.^2-GausCurv); % larger principal curvature
            
            
            %% Report time
            measurement.k1 = k1;
            measurement.k2 = k2;
            measurement.H_rec = MeanCurv;
            measurement.H_rangeInPatch = H_rangeInPatch;
            measurement.K_rec = GausCurv;
            measurement.fittingParameters = fittingParameters;
            measurement.a_patch = a_patch;
            measurement.b_patch = b_patch;
            measurement.coordsRecSurface = coordsRecSurface;
            measurement.surfaceNormalsPoly3 = surfaceNormalsPoly3;
            measurement.MeanCurvStdErr = MeanCurvStdErr;
            measurement.ChiSqrd = ChiSqrd;
            measurement.Residuals = Residuals;
            measurement.LengthScaleVariation = LengthScaleVariation;
            measurement.LengthScalePatch = LengthScalePatch;
            iterNum = measurement.iterNum;
            if iterNum == 1
                %circular patch fits
                measurement.Alphas = nan(measurement.numberOfPoints,...
                    measurement.numOfIterations);
                measurement.Betas = nan(measurement.numberOfPoints,...
                    measurement.numOfIterations);
                measurement.relErrH = nan(measurement.numberOfPoints,...
                    measurement.numOfIterations);
            end
            measurement.Alphas(:,iterNum) = Alphas;
            measurement.Betas(:,iterNum) = Betas;
            %relErrH = (MeanCurv - measurement.H_anal)./measurement.H_anal;
            %measurement.relErrH(:,iterNum) = relErrH;
            measurement.GOF.h = GOF_h;
            measurement.GOF.p = GOF_p;
            measurement.P11F = P11F;
            measurement.P20F = P20F;
            measurement.P02F = P02F;
            measurement.NumPerPatch = NumPerPatch;
            
        end
        function [MeanCurv,GausCurv,patchRadii,...
                ResidualsInnerMean,ResidualsOuterMean,...
                ResidualsTimesRadius,H_rangeInPatch,...
                fittingParameters,a_patch,b_patch,...
                coordsRecSurface,...
                surfaceNormalsPoly3,MeanCurvStdErr,...
                ChiSqrd,LengthScaleVariation,...
                LengthScalePatch,Alphas,Betas,GOF_h,...
                GOF_p,P11F,P20F,P02F,NumPerPatch] = ...
                littleAnalysisLoop(i,CoordsIn,rp,NS_kd,k_est,minRp,maxHVar,center,d_lambda)
            %% GET CIRCULAR PATCH COORDINATES
            queryPoint = CoordsIn(i,:);
            patchTooBig = 1;
            while patchTooBig
                [idx,~]=rangesearch(NS_kd,queryPoint,rp); % this delivers the indices of points in the neighborhood
                patchRadii = rp;
                idx = idx{1};
                idx=nonzeros(idx);
                %this is the neighborhood
                patchCoordsGF=CoordsIn(idx,:);
                usefulPoints = length(patchCoordsGF(:,1));
                H_inPatch = k_est(idx);
                deltaH_rec = range(H_inPatch);
                if (rp==minRp)||(usefulPoints==6)
                    patchTooBig = 0;
                    if usefulPoints<6
                        [idx,dst]=knnsearch(NS_kd,queryPoint,'k',6); % this delivers the indices of points in the neighborhood
                        idx=nonzeros(idx);
                        %this is the neighborhood
                        patchCoordsGF=CoordsIn(idx,:);
                        patchRadii = max(dst(:));
                        usefulPoints = length(patchCoordsGF(:,1));
                    end
                else
                    patchTooBig = (deltaH_rec/k_est(i))>maxHVar;
                    if patchTooBig
                        rp = rp*(1-2*pi/d_lambda);
                        if (rp<minRp)||(usefulPoints<6)
                            rp=minRp;
                        end
                    end
                end
            end
            if (usefulPoints<6)
                error('Too few useful points.')
            end
            [eigenvectors,~] = DropRecon3Dv1.getPrincipalAxes(...
                queryPoint,patchCoordsGF);
            % e1 cooresponds to normal
            e1=eigenvectors(:,1);e2=eigenvectors(:,2);...
                e3=eigenvectors(:,3);
            % fix orientation of normal
            if sign((center - queryPoint)*e1) == -1
                e1 = -e1;
                e2 = -e2;
                %e3 = e3;
            end
            
            if (sum([isnan(e1)', isnan(e2)', isnan(e3)'])>0)
                warning('Principal axes are meaningless')
            end
            thisCoordinate = [nan, nan, nan];
            gofH = nan;
            gofP = nan;
            p20f = nan;
            p02f = nan;
            p11f = nan;
            N_GF = [nan, nan, nan];
            usefulPoints = length(patchCoordsGF(:,1));
            if (usefulPoints<6)
                warning('Too few points for meaningful fit')
                H = nan;
                dH = nan;
                K = nan;
                inRes = nan;
                outRes = nan;
                RRes = nan;
                LengthVar = nan;
                LengthPatch = nan;
                fPar2 = nan(1,6);
                chiSq2 = nan;
                alpha = nan;
                beta = nan;
                npp = usefulPoints;
            else
                [patchCoordsLF,W_patch,qPoint] ...
                    = DropRecon3Dv1.getCoordsLocalFrame(queryPoint,patchCoordsGF,e1,e2,e3 );
                %% FIT CIRCULAR PATCH COORDINATES
                X_patch=patchCoordsLF(:,1);
                Y_patch=patchCoordsLF(:,2);
                Z_patch=patchCoordsLF(:,3);
                
                [fPar2,errfPar2,chiSq2,residuals2] = ...
                    DropRecon3Dv1.getPolyFit2(X_patch,...
                    Y_patch,Z_patch,W_patch);
                
                p20f = fPar2(1);
                p02f = fPar2(2);
                p11f = fPar2(3);
                p10 = fPar2(4);
                p01 = fPar2(5);
                p00 = fPar2(6);
                
                qx0=qPoint(1);qy0=qPoint(2);
                % This is how we turn the fit data into mean curvature
                
                [H,dH]=DropRecon3Dv1.getMeanCurvaturePoly2(...
                    fPar2,qx0,qy0,errfPar2);
                % This is how we turn the fit data into gaussian curvature
                [K]=DropRecon3Dv1.getGaussianCurvaturePoly2(...
                    fPar2,qx0,qy0);
                
                [nxLF,nyLF,nzLF]=DropRecon3Dv1.getNormalPoly2(...
                    fPar2,qx0,qy0);
                N_LF = [nxLF,nyLF,nzLF];
                
                alpha = std(residuals2);
                beta = mean(residuals2);
                [gofH,gofP] = chi2gof(residuals2);
                R_patch = sqrt(X_patch.^2+Y_patch.^2+...
                    Z_patch.^2);
                R_median = median(R_patch);
                inRes = mean(residuals2(R_patch<R_median));
                outRes = mean(residuals2(R_patch>R_median));
                RRes = sum(residuals2.*R_patch)./sum(R_patch);
                LengthVar = nan;
                LengthPatch = patchRadii;
                npp = length(X_patch);
                qzFit = p20f*qx0^2+p02f*qy0^2+...
                    p11f*qx0*qy0+p10*qx0+p01*qy0+p00;
                qPointLF = [qx0,qy0,qzFit];
                N_GF = N_LF/[e2,e3,e1];
                thisCoordinate = qPointLF/[e2,e3,e1]+...
                    mean(patchCoordsGF);
            end
            coordsRecSurface = thisCoordinate;
            surfaceNormalsPoly3 = N_GF;
            H_rangeInPatch = deltaH_rec;
            MeanCurv = H;
            MeanCurvStdErr = dH;
            GausCurv = K;
            ChiSqrd = chiSq2;
            ResidualsInnerMean = inRes;
            ResidualsOuterMean = outRes;
            ResidualsTimesRadius = RRes;
            LengthScaleVariation = LengthVar;
            LengthScalePatch = LengthPatch;
            a_patch = LengthPatch;
            b_patch = LengthPatch;
            fittingParameters = fPar2;
            Alphas = alpha;
            Betas = beta;
            GOF_h = gofH;
            GOF_p = gofP;
            P11F = p11f;
            P02F = p02f;
            P20F = p20f;
            NumPerPatch = npp;
            
        end
        function [CoordsLocalFrame,distWeights,qPointLF] = ...
                getCoordsLocalFrame(qPoint,nnPoints,e1,e2,e3)
            gaussianWeighted=0; % Flip to 1 if you want coordinates gaussian weighted
            %%
            %   project nnPoints onto tangent plane
            %   defined by the eigenvectors
            %%
            
            nofnnPoints=length(nnPoints(:,1));
            
            %subtract query point from all neighborhood points
            nnPointsqPoint=nnPoints-ones(nofnnPoints,1)*qPoint;
            
            %calculate distances from neighborhood points to query point
            d=sqrt(nnPointsqPoint(:,1).^2+nnPointsqPoint(:,2).^2+nnPointsqPoint(:,3).^2);
            
            %create length scale for gaussian weighting. h will be ? standard
            %deviations
            h=max(d);
            %center of mass
            COM=mean(nnPoints);
            
            %project query point onto tangent plane which goes through center of mass
            qx0=(qPoint-COM)*e2;
            qy0=(qPoint-COM)*e3;
            qz0=0;
            qPointLF=[qx0,qy0,qz0];
            
            
            %coordinates with com subtracted
            nnPointsM = nnPoints-ones(nofnnPoints,1)*COM;
            % Matrix to rotate coords into "Local Frame"
            nn2LF_Mat = [e2,e3,e1];
            %coordinates in local frame
            CoordsLF = nnPointsM*nn2LF_Mat;
            %calculate weights for points in neighborhood
            
            
            if gaussianWeighted
                CoordsLFqPointLF=CoordsLF-ones(nofnnPoints,1)*qPointLF;
                dsq=sqrt(CoordsLFqPointLF(:,1).^2+CoordsLFqPointLF(:,2).^2+CoordsLFqPointLF(:,3).^2);
                distWeights=(exp(-1/2*dsq.^2./h^2))';
            else
                distWeights=ones(nofnnPoints,1);
            end
            %scatter3(CoordsLF(:,1),CoordsLF(:,2),CoordsLF(:,3),1,distWeights);
            CoordsLocalFrame=CoordsLF;
            %close all
            
        end
        function [fPar2,errfPar2,chsq2,residuals2] = getPolyFit2(X,Y,Z,W)
            XFIT=X;
            YFIT=Y;
            ZFIT=Z;
            WFIT1=ones(length(XFIT),1);
            WFIT2=W;
            
            %% 0. fit plane
            
            mx=(max(ZFIT)-min(ZFIT) )/(max(XFIT)-min(XFIT) );
            my=(max(ZFIT)-min(ZFIT) )/(max(YFIT)-min(YFIT) );
            p10i=mx;p01i=my; p00i=mean(ZFIT);p00i=0;
            
            px=min(XFIT):(max(XFIT)-min(XFIT))/20:max(XFIT);
            py=min(YFIT):(max(YFIT)-min(YFIT))/20:max(YFIT);
            
            [X0,Y0]=meshgrid(px,py);
            x0=X0(:);y0=Y0(:);
            
            startP0=[p10i,p01i,p00i];
            
            %% 1. fit deg 1 polynomial in X,Y
            
            [fParOut0,~,~]=fitPoly1P10P01P00(startP0,XFIT,YFIT,ZFIT,WFIT1);
            p10f=fParOut0(1);p01f=fParOut0(2);p00f=fParOut0(3);
            
            z0=p10f*x0+p01f*y0+p00f;
            Z0=reshape(z0,size(X0));
            
            %% 2. fit deg 2 polynomial in X
            p11i=0 ;p10i= p10f;p01i= p01f;p00i= p00f;
            
            p20i=0;p02i=0;
            startP1=[p20i,p11i,p10f,p01f,p00f];
            [fParOut1,~,~]=fitPoly2X(startP1,XFIT,YFIT,ZFIT,WFIT2);
            
            %% 3. fit deg 2 polynomial in X,Y
            
            startP2=[p20i fParOut1];
            [fParOut2,errfParOut2,chsq2Out]=fitPoly2(startP2,XFIT,YFIT,ZFIT,WFIT2);
            Z_fitted = @(x,y) fParOut2(1)*x.^2 + fParOut2(2)*y.^2 + ...
                fParOut2(3)*x.*y + fParOut2(4)*x + fParOut2(5)*y + ...
                fParOut2(6);
            
            fPar2=fParOut2;
            errfPar2=errfParOut2;
            chsq2=chsq2Out;
            residuals2 = Z_fitted(XFIT(:),YFIT(:)) - ZFIT(:);
            
            
            needed = 0;
            if needed
                Z2 = Z_fitted(X0(:),Y0(:));
                Z2=reshape(Z2,size(X0));
                figure(1)
                hold off;
                scatter3(XFIT,YFIT,ZFIT)
                hold on;
                mesh(X0,Y0,Z2)
                hold off;
            end
        end
        function [H,dH] =getMeanCurvaturePoly2(fpar2,x,y,errfpar2)
            %% FIT PARAMETERS, assigned
            p30 = 0;
            p03 = 0;
            p21 = 0;
            p12 = 0;
            p20 = fpar2(1);
            p02 = fpar2(2);
            p11 = fpar2(3);
            p10 = fpar2(4);
            p01 = fpar2(5);
            %p00 = fpar(6);
            
            %% Derivatives of height with respect to x,y (parameterized as u,v, for visual clarity)
            hu = p10 + 2*p20*x + 3*p30*x^2 + p11*y + 2*p21*x*y + p12*y^2;
            hv = p01 + p11*x + p21*x^2 + 2*p02*y + 2*p12*x*y + 3*p03*y^2;
            huu = 2*p20 + 6*p30*x + 2*p21*y;
            hvv = 2*p02 + 2*p12*x + 6*p03*y;
            huv = p11 + 2*p21*x + 2*p12*y;
            
            %% Mean Curvature, H
            H = ((1+hv^2)*huu-2*hu*hv*huv+(1+hu^2)*hvv)./(2*(1+hu^2+hv^2)^(3/2));
            if nargin>3
                % Error Propogation
                % Standard Errors of Fit Parameters
                dp30 = 0;
                dp03 = 0;
                dp21 = 0;
                dp12 = 0;
                dp20 = errfpar2(1);
                dp02 = errfpar2(2);
                dp11 = errfpar2(3);
                dp10 = errfpar2(4);
                dp01 = errfpar2(5);
                %dp00 = errfpar(6);
                % Derivatives of hu, hv, huu, & hvv w.r.t fit parameters.
                % hu
                hu_p10 = 1;
                hu_p20 = 2*x;
                hu_p30 = 3*x^2;
                hu_p11 = y;
                hu_p21 = 2*x*y;
                hu_p12 = y^2;
                % hv
                hv_p01 = 1;
                hv_p11 = x;
                hv_p21 = x^2;
                hv_p02 = 2*y;
                hv_p12 = 2*x*y;
                hv_p03 = 3*y^2;
                % huu
                huu_p20 = 2;
                huu_p30 = 6*x;
                huu_p21 = 2*y;
                % hvv
                hvv_p02 = 2;
                hvv_p12 = 2*x;
                hvv_p03 = 6*y;
                % huv
                huv_p11 = 1;
                huv_p21 = 2*x;
                huv_p12 = 2*y;
                
                % Standard Errors of hu, hv, huu, & hvv:
                dhu = sqrt((hu_p10*dp10)^2+(hu_p20*dp20)^2+(hu_p30*dp30)^2+...
                    (hu_p11*dp11)^2+(hu_p21*dp21)^2+(hu_p12*dp12)^2);
                dhv = sqrt((hv_p01*dp01)^2+(hv_p11*dp11)^2+(hv_p21*dp21)^2+...
                    (hv_p02*dp02)^2+(hv_p12*dp12)^2+(hv_p03*dp03)^2);
                dhuu = sqrt((huu_p20*dp20)^2+(huu_p30*dp30)^2+...
                    (huu_p21*dp21)^2);
                dhvv = sqrt((hvv_p02*dp02)^2+(hvv_p12*dp12)^2+...
                    (hvv_p03*dp03)^2);
                dhuv = sqrt((huv_p11*dp11)^2+(huv_p21*dp21)^2+...
                    (huv_p12*dp12)^2);
                % Derivatives of H w.r.t. hu, hv, huu, & hvv
                H_hu = (-2*huv*hv + 2*hu*hvv)/(2*(1 + hu^2 + hv^2)^(3/2)) -...
                    ( 3*hu*(-2*hu*huv*hv + huu*(1 + hv^2) + (1 + hu^2)*hvv))/...
                    ( 2*(1 + hu^2 + hv^2)^(5/2));
                H_hv = (-2*hu*huv + 2*huu*hv)/(2*(1 + hu^2 + hv^2)^(3/2)) -...
                    ( 3*hv*(-2*hu*huv*hv + huu*(1 + hv^2) + (1 + hu^2)*hvv))/...
                    ( 2*(1 + hu^2 + hv^2)^(5/2));
                H_huu = (1 + hv^2)/(2*(1 + hu^2 + hv^2)^(3/2));
                H_hvv = (1 + hu^2)/(2*(1 + hu^2 + hv^2)^(3/2));
                H_huv = -((hu*hv)/(1 + hu^2 + hv^2)^(3/2));
                % Standard Error of H
                dH=sqrt((H_hu*dhu)^2+(H_hv*dhv)^2+(H_huu*dhuu)^2+...
                    (H_hvv*dhvv)^2+(H_huv*dhuv)^2);
            else
                dh=nan;
            end
            
        end
        function [K] =getGaussianCurvaturePoly2(fpar2,x,y)
            %% FIT PARAMETERS, assigned
            p30 = 0;
            p03 = 0;
            p21 = 0;
            p12 = 0;
            p20 = fpar2(1);
            p02 = fpar2(2);
            p11 = fpar2(3);
            p10 = fpar2(4);
            p01 = fpar2(5);
            %p00 = fpar(6);
            
            %% Derivatives of height with respect to x,y (parameterized as u,v, for visual clarity)
            hu = p10 + 2*p20*x + 3*p30*x^2 + p11*y + 2*p21*x*y + p12*y^2;
            hv = p01 + p11*x + p21*x^2 + 2*p02*y + 2*p12*x*y + 3*p03*y^2;
            huu = 2*p20 + 6*p30*x + 2*p21*y;
            hvv = 2*p02 + 2*p12*x + 6*p03*y;
            huv = p11 + 2*p21*x + 2*p12*y;
            
            %% Gaussian Curvature, K
            K = (huu*hvv - huv^2)./((1+hu^2+hv^2)^2);
            
        end
        function [nxLF,nyLF,nzLF] = getNormalPoly2(fpar2,xq,yq)
            %% FIT PARAMETERS, assigned
            p30 = 0;
            p03 = 0;
            p21 = 0;
            p12 = 0;
            p20 = fpar2(1);
            p02 = fpar2(2);
            p11 = fpar2(3);
            p10 = fpar2(4);
            p01 = fpar2(5);
            %% Derivatives of function with respect to x,y,
            % evaluated at xq, yq
            fx = p10 + 2*p20*xq + 3*p30*xq^2 + p11*yq + 2*p21*xq*yq + p12*yq^2;
            fy = p01 + p11*xq + p21*xq^2 + 2*p02*yq + 2*p12*xq*yq + 3*p03*yq^2;
            
            N = [fx,fy,-1];
            nxLF = N(1)/norm(N);
            nyLF = N(2)/norm(N);
            nzLF = N(3)/norm(N);
            
        end
        function [measurement] = applyMedianFilter2Measurement(measurement)
            %
            if isempty(measurement.medFiltPatch.idx)
                CoordsIn = [measurement.x,measurement.y,measurement.z];
                NS_kd = KDTreeSearcher(CoordsIn);
                startIndex = 1;
                endIndex = length(measurement.x);
                patchSize = measurement.medianFilterSizeInPixels;
                H_in = measurement.H_rec;
                H_out = nan(endIndex,1);
                for i = startIndex:endIndex
                    %% GET CIRCULAR PATCH COORDINATES
                    queryPoint = CoordsIn(i,:);
                    
                    [idx,~]=rangesearch(NS_kd,queryPoint,patchSize); % this delivers the indices of points in the neighborhood
                    idx = idx{1};
                    idx=nonzeros(idx);
                    measurement.medFltrPatch.idx{i} = idx;
                    %this is the neighborhood
                    H_patch = H_in(idx);
                    H_out(i) = median(H_patch);
                end
            else
                startIndex = 1;
                endIndex = length(measurement.x);
                H_in = measurement.H_rec;
                H_out = nan(endIndex,1);
                for i = startIndex:endIndex
                    %% GET CIRCULAR PATCH COORDINATES
                    idx = measurement.medFiltPatch.idx{i};
                    %this is the neighborhood
                    H_patch = H_in(idx);
                    H_out(i) = median(H_patch);
                end
            end
            measurement.H_rec = H_out;
            [measurement] = DropRecon3Dv1.applyMedianFilter2Recon(measurement);
            %
        end
        function [measurement] = findMedianHinPatch(measurement)
            CoordsIn = [measurement.x,measurement.y,measurement.z];
            NS_kd = KDTreeSearcher(CoordsIn);
            startIndex = 1;
            endIndex = length(measurement.x);
            patchSize = measurement.a_patch;
            H_rec = measurement.H_rec;
            H_anal = measurement.H_anal;
            H_medRecInPatch = nan(endIndex,1);
            H_medAnalInPatch = nan(endIndex,1);
            for i = startIndex:endIndex
                %% GET CIRCULAR PATCH COORDINATES
                queryPoint = CoordsIn(i,:);
                
                [idx,~]=rangesearch(NS_kd,queryPoint,patchSize(i)); % this delivers the indices of points in the neighborhood
                idx = idx{1};
                idx=nonzeros(idx);
                measurement.medFltrPatch.idx{i} = idx;
                %this is the neighborhood
                H_recPatch = H_rec(idx);
                H_medRecInPatch(i) = median(H_recPatch);
                H_analPatch = H_anal(idx);
                H_medAnalInPatch(i) = median(H_analPatch);
            end
            measurement.H_medRecInPatch = H_medRecInPatch;
            measurement.H_medAnalInPatch = H_medAnalInPatch;
        end
        function [] = makeGif(X,Y,Z,V,cAxisStr,gifPath)
            % 30 July 2017
            numOfFrames = 192;
            close all
            
            colormap('parula')
            
            hSurface = surf(X,Y,Z,V,'LineStyle','none');
            axis equal
            view([0 -180]);
            %colorbar;
            
            set(gca,'color',[0 0 0])
            set(gcf,'color',[0 0 0])
            %title(cAxisStr,'Color','white')
            set(gcf, 'InvertHardCopy', 'off');
            set(gca,'nextplot','replacechildren','visible','off')
            f = getframe;
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,numOfFrames) = 0;
            X= [zeros(numOfFrames/4,1);ones(numOfFrames/4,1);...
                zeros(numOfFrames/4,1);ones(numOfFrames/4,1)];
            Y = zeros(numOfFrames,1);
            Z = [ones(numOfFrames/4,1);zeros(numOfFrames/4,1);...
                ones(numOfFrames/4,1);zeros(numOfFrames/4,1)];
            direction = [X Y Z];
            for k = 1:numOfFrames
                rotate(hSurface,direction(k,:),2*(360/numOfFrames))
                f = getframe;
                A = rgb2ind(f.cdata,map,'nodither');
                [m,n] = size(A);
                im(1:m,1:n,1,k) = A;
            end
            imwrite(im,map,gifPath,'DelayTime',0,...
                'LoopCount',inf,'WriteMode','overwrite')
            
            
            % SAVE COLORMAP
            
            %f=figure;
            %a=axes;
            set(gca,'color',[1 1 1])
            set(gcf,'color',[1 1 1])
            c=colorbar;
            ylabel(c,cAxisStr) 
            epsPath = [gifPath(1:end-4),'-colormap.eps'];
            print(epsPath,'-depsc')
            
            close all
            clear im
            
            
        end
        
        function [evalMetric,measurement] = evalMeas(measurement,anal,evalMetric)
            %% Make Gifs and Maps
            measurement.H_analGrid = anal.HGrid;
            measurement.H_anal = anal.H(:);
            measurement.s_anal = anal.s(:);
            measurement.relErrH = ...
                (measurement.H_rec-measurement.H_anal)./measurement.H_anal;
            measurement.relErrHStar = ...
                (measurement.H_rec-measurement.H_anal)./measurement.H_rec;
            
            evalMetric.HRelErr = measurement.relErrH;
            evalMetric.HErr = (measurement.H_rec-measurement.H_anal);
            
            evalMetric.HRelErr_mean = mean(evalMetric.HRelErr);
            evalMetric.HErr_mean = mean(evalMetric.HErr);
            
            evalMetric.HRelErr_std = std(evalMetric.HRelErr);
            evalMetric.HErr_std = std(evalMetric.HErr);
            
            measurement = DropRecon3Dv1.makeAnimatedGifs(measurement);
            DropRecon3Dv1.makeProjectionMap(measurement)
            
            
            hold off
            figure(2)
            hist(measurement.relErrH);
            title('Histogram of Relative Errors of H')
            xlabel('Relative Error of Mean Curvature')
            ylabel('Counts')
            
            fname = [measurement.path2ResultsDir,'RelErrH.eps'];
            print(gcf,fname,'-depsc');
        end
    end
    
end


