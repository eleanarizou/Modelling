function [RealradialAvgNuc] = ModRealColoniesAnalysis(meta,dataDir,fInFormat, MAINFILENAME, thres, FILLs,chans,DAPIChannel)
RealradialAvgNuc = []
for jjj = 1:length(FILLs)
    
    meta.channelLabel = {chans{jjj}{1},chans{jjj}{2},'DAPI'};
    save(fullfile(dataDir,'metaData.mat'),'meta');
    
    %% process (save DAPI channel MIP for segmentation)
    
    % although findColonies can find multiple colonies in a single image,
    % for the LSM there is one per image, which this code below assumes
    findColParam = struct('sclose',10, 'sopen', 4, 'checkcontained', false,...
        'minArea', [],'convhull', true);
    close all;
    colonies(numel([FILLs{jjj}{:}])) = Colony;
    DAPnorm = {true, false};
    DAPnormFilname = {['_'],['noDAPI']};
    for u = 1:2
        for coli = [FILLs{jjj}{:}]
            
            % cleanScale in micron
            param = {'DAPIChannel',DAPIChannel, 'colID',coli, 'adjustmentFactor', thres{jjj},...
                'prenormImage',DAPnorm{u},'clparameters',findColParam};
            
            filename = sprintf(fInFormat, coli);
            colony = ModprocessOneColonyImage(filename, dataDir, param);
            colonies(coli) = colony;
            colonies(coli).setID(coli);
            % setID overwrites filename, which was the right thing for epi, not
            % here
            colonies(coli).filename = filename;
        end
        save(fullfile('/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling','colonies'), 'colonies');
        
        %% show plot of different conditions side by side
        load(fullfile('/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling','colonies'), 'colonies');
        nchan = 2;
        set(gcf,'Position',[0 0 800 1000])
        aveTim = []
        for ii = 1:length(FILLs{jjj})
            colsNow = colonies(FILLs{jjj}{ii});
            avgstruct = makeAveragesNoSegmentation(meta,350,3,colsNow);
            for jj = 1:nchan
                q = (ii-1)*nchan+jj;
                subplot(length(FILLs{jjj}),nchan,q);
                plot(avgstruct.r,avgstruct.nucAvg(:,jj),'k.-','LineWidth',3); hold on;
                %                 ylim([0 350])
                for kk = 1:length(colsNow)
                    %all the colonies radial profile and the average trend
                    plot(avgstruct.r,colsNow(kk).radialProfile.NucAvg(:,jj),'c-*','LineWidth',3);
                    %                     ylim([0 350])
                end
                title([meta.conditions{ii}, ' Files ' ,num2str(FILLs{jjj}{ii}),' ', chans{jjj}{ii,jj}]); hold off;
            end
            aveTim = [aveTim,{avgstruct}]
        end
        set(gcf,'Position',[0 0 1100 1500])
        saveas(gcf,fullfile('/Volumes/storage/Eleana/modelling_gastruloids/XMASmodellling',[MAINFILENAME, DAPnormFilname{u},chans{jjj}{ii,1},'.png']));
        close all;
        
    end
        RealradialAvgNuc = [RealradialAvgNuc,{aveTim}]
end
    
end
