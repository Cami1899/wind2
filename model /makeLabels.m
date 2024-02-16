%Function for making plot labels easily
function [] = makeLabels(thisTitle, thisxlabel, thisylabel)
        title(thisTitle, 'FontSize', 20)
        xlabel(thisxlabel, "FontSize", 15)
        ylabel(thisylabel, "FontSize", 15)
        set(gcf,'Position',[0 0 1500 500])
end