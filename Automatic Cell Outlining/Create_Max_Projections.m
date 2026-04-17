folderPath = '\\10.67.15.232\lab\Morgane\MV635-...-642_26-03-12_IF NBDY Dox+ - _ Dcpa1a 647 LSM14A 488';
identifiers = {'470S','395S','555S','640S'};

maxProjectTiffFolder(folderPath, identifiers);

function maxProjectTiffFolder(folderPath, identifiers)

    warning('off', 'MATLAB:imagesci:tiff:libraryWarning');
    folderPath  : 'Z:\Marius\20260327_WT-70ko_SunTag_IF_smFISH_Hip-Tor\20260327'
    identifiers : {'555S','395S', '640S'}
    
    % Create output folder
    outputFolder = fullfile(folderPath, 'Max Projections');
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    
    % Get all tif and tiff files
    files = dir(fullfile(folderPath, '*.tif'));
    files = [files; dir(fullfile(folderPath, '*.tiff'))];
    
    if isempty(files)
        error('No TIFF files found in folder.');
    end
    
    % Loop through files
    for i = 1:numel(files)
    
        fname = files(i).name;
    
        % Check for identifiers
        if ~any(contains(fname, identifiers))
            continue
        end
    
        inputPath = fullfile(folderPath, fname);
    
        % Output filename
        [~, baseName, ~] = fileparts(fname);
        outputPath = fullfile(outputFolder, [baseName '_MAX.tif']);
    
        % ---- Read Z-stack efficiently ----
        t = Tiff(inputPath, 'r');
    
        % Read first slice
        maxProj = t.read();
    
        % Loop through Z slices
        while ~t.lastDirectory()
            t.nextDirectory();
            slice = t.read();
            maxProj = max(maxProj, slice);
        end
    
        t.close();
    
        % ---- Save max projection ----
        imwrite(maxProj, outputPath, 'Compression', 'none');
    
        fprintf('Saved max projection: %s\n', outputPath);
    end
    
    fprintf('All matching files processed and saved to z-stacks.\n');

end

