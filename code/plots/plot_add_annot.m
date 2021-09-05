function [ax,legendFig,colorAll] = plot_add_annot(axisP,annotStruct,isRow,legendAxis,inOpts)


    defaultOpts.legendFig = 1;
    defaultOpts.sortLegend = 1;
    defaultOpts.legendOrient = 'Vertical';


    colorSet = @(zCol)cell2mat(cellfun(@(x)hexColor(x(2:end)),zCol,'uniformoutput',0))';
    % defaultOpts.groupColorMap = colorSet({'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#117744', '#44AA77', '#88CCAA', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788', '#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788' });

    defaultOpts.groupColorMap{1} = brewermap(12,'Dark2');
    defaultOpts.groupColorMap{2} = brewermap(8,'Set1');
    defaultOpts.groupColorMap{3} = brewermap(12,'Set3');

    defaultOpts.grey = [0.6627 0.6627 0.6627];

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts)

    % ax = subplot('Position',axisPos);

    cfnames = fieldnames(annotStruct);
    colIdx = 0;
    for i = 1:length(cfnames)

        cf = cfnames{i};

        annotB = annotStruct.(cfnames{i});
        [groupIDidx,gName] = grp2idx(annotB);

        if opts.sortLegend
            if isnumeric(annotB)
                [gName,zidxFix] = sort(str2double(gName));
            else
                [gName,zidxFix] = sort(gName);
            end
            [~,zidxFix] = sort(zidxFix);

            groupIDidx  = zidxFix(groupIDidx) + colIdx;
        end

        colIdx = colIdx + length(gName);

        if isRow
            groupID{i,1} = groupIDidx(:)';
        else
            groupID{1,i} = groupIDidx(:);
        end

        if isfield(opts,[ cf '_cmap' ])
            cmap = opts.([ cf '_cmap' ]);
        else
            if iscell(opts.groupColorMap)
                cmap = opts.groupColorMap{mod(i-1,length(opts.groupColorMap))+1};
            else
                cmap = opts.groupColorMap;
            end
        end

        if isa(cmap,'containers.Map')
            cfCmap = nanvalues(cmap,gName);

            cfCmap(isemptycell(cfCmap)) = { opts.grey };
            cfCmap = cell2mat(cfCmap);
        else
            cfCmap = cmap(mod((1:length(gName))-1,length(cmap))+1,:);
        end

        gNameAll{i,1} = gName(:);
        colorAll{i,1} = cfCmap;
    end


    groupLen = cellfun(@(x)length(x),gNameAll);
    groupID = cell2mat(groupID);
    colorM = cell2mat(colorAll);

    iax = image(axisP,groupID);
    ax = iax.Parent;
    
    colormap(ax,colorM);

    if isRow
        ax.YTick = 1:length(cfnames);
        ax.YTickLabels = cfnames;
        ax.XTick = [];        
    else
        ax.XTick = [];        
        ax.YTick = [];
    end
    
    legendFig = [];
    if opts.legendFig
        
        % Configures axis
        if ~exist('legendAxis','var') || isempty(legendAxis)
            legendFig = figure();
            
            lcf = length(cfnames);
            if lcf == 1
                legendAxis{1} = gca();
            elseif lcf < 3
                legendAxis = arrayfun(@(x)subplot(1,lcf,x),1:lcf,'unif',0);
            else
                legendAxis = arrayfun(@(x)subplot(2,ceil(lcf/2),x),1:lcf,'unif',0);
            end

        elseif iscell(legendAxis)            
            %lcf = length(legendAxis);
            lcf = min(length(legendAxis),length(groupLen));
            if (lcf < length(cfnames))
                warning('Insufficient number of axis to plot all annotations');
            end                          
        elseif isa(legendAxis,'matlab.graphics.axis.Axes')
            lcf = min(length(legendAxis),length(groupLen));
            legendAxis = { legendAxis };
        end
        
        cidx = 1;
        qq = cell(lcf,1);
        leg = cell(lcf,1);      

        for zi = 1:lcf
            cAx = legendAxis{zi};
            axes(cAx);
            zp = cidx:(cidx+groupLen(zi)-1);
            
            for ii = 1:length(zp)
                qq{zi}(ii) = patch(NaN, NaN, colorM(zp(ii),:));
            end

            cidx = cidx + groupLen(zi);

            % axL = gca;
            gName = gNameAll{zi};
            gName = gName(:);
            
            if isnumeric(gName)
                leg{zi} = legend(qq{zi},num2cellstr(gName),'Location','NorthWest','Orientation',opts.legendOrient);
            else
                leg{zi} = legend(qq{zi},gName,'Location','NorthWest','Orientation',opts.legendOrient);
                leg{zi}.Interpreter = 'None';                
            end
            
            cAx.Visible = 'Off';
            disp('Ax');
        end
    end



end

%         iax = image(axisP,groupIDidx);
%         ax = iax.Parent;
%
%         ax.YTick = [];
%         ax.XTick = [];
%
%
%         colormap(ax,cmap);
%
%         if opts.legendFig
%             legendFig = figure();
%
%             for ii = 1:size(cmap,1)
%                 qq(ii) = patch(NaN, NaN, cmap(ii,:));
%             end
%
%             axL = gca;
%             if isnumeric(gName)
%                 leg = legend(qq,num2cellstr(gName),'Location','NorthWest','Orientation','horizontal');
%             else
%                 leg2 = legend(qq,gName,'Location','NorthWest','Orientation','horizontal');
%                 leg2.Interpreter = 'None';
%             end
%             axL.Visible(1);
%         end
%     end
