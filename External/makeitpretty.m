if exist('magnification')==0
    magnification=1;
end
ff=gcf;
ff.Color='White';
for i=1:numel(ff.Children)
    ff.Children(i).FontSize=14*magnification;
    ff.Children(i).FontName='Roman';
    try
    ff.Children(i).Box='on';
    end
    try
        ff.Children(i).Interpreter='latex' ;
    end
    try
        ff.Children(i).TickLabelInterpreter='latex' ;
    end
    try
        ff.Children(i).XLabel.Interpreter='latex';
        ff.Children(i).XLabel.FontSize=18*magnification;
        ff.Children(i).YLabel.Interpreter='latex';
        ff.Children(i).YLabel.FontSize=18*magnification;
        try
            ff.Children(i).ZLabel.Interpreter='latex';
            ff.Children(i).ZLabel.FontSize=18*magnification;
        end
    end
    try   
        ff.Children(i).Title.Interpreter='latex';
        ff.Children(i).Title.FontSize=20*magnification;
    end
    try
        ff.CurrentAxes=ff.Children(i); axis tight 
        ff.Children(i).YLim=ff.Children(i).YLim+[-0.05 0.05]*(ff.Children(i).YLim(2)-ff.Children(i).YLim(1)); %increases y axes by 1%
        ff.Children(i).XLim=ff.Children(i).XLim+[-0.01 0.01]*(ff.Children(i).XLim(2)-ff.Children(i).XLim(1)); %increases y axes by 1%
        try
                    ff.Children(i).ZLim=ff.Children(i).ZLim+[-0.01 0.01]*(ff.Children(i).ZLim(2)-ff.Children(i).ZLim(1)); %increases y axes by 1%

        end
    end
    
    
%     for k=1:numel(ff.Children(i).Children)
%         try
%             ff.Children(i).Children(k).LineWidth=1.2*magnification;
%         end
%     end
%         
        
end

