function [M,KM, SK] = Nutrient(path,titles)
 arguments 
        path="datasets\Nutrient.mat"; 
        titles = ["Calcium", "Fer", "Proteine","Vitamine A", "Vitamine C"];
end 
    if isfile(path)
       Data= load(path);
    else
    disp([path, "Not Found"]);
    end 
    M = VectorNormalised(VectorOfRanks(Data.Nutrient));
    fig = figure;
    [~,ax,BigAx,H,HAx] = plotmatrix(M);
    savefig("outputs/figures/scatterplotnutrient.fig");
    saveas(fig,"outputs/figures/scatterplotnutrient.png","png");
    KM = Matrix_KendallSpearman(M);
    iterations = size(ax,1);
    SK = skewness(VectorOfRanks(Data.Nutrient));
   
    for i = 1:iterations
        
        ax(i,1).YLabel.String = titles(i);
        ax(iterations,i).XLabel.String = titles(i);
    end
    
    %close;
end