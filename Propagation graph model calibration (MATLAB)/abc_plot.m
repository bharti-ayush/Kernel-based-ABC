% Function to plot the density estimates of the posterior
function [] = abc_plot(samples, xlabels, xlimits, true_param, fig_number)
    if length(samples(1,:)) == 5
        nrow = 3;
    else
        nrow = 2;
    end
    
    ncol = 2;
    figure(fig_number);
%     for ii=1:length(xlabels)
    for ii=1:nrow.*ncol
%         [f,xi] = ksdensity(samples(:,ii));
%         subplot(1,length(xlabels),ii)
        hsub = subplot(nrow, ncol, ii);
        if ii == 2
            [f,xi] = ksdensity(samples(:,ii));
           histogram( round(samples(:,ii) ) , 'FaceColor', 'black', "Normalization", "probability", "EdgeColor", "black");
           xlim(xlimits(ii,:));
           line([round(mean(samples(:,ii))), round(mean(samples(:,ii)))], [0, max(f)] , 'Color', 'r', 'linewidth', 1.5);
           line([true_param(ii), true_param(ii)], [0, max(f)] , 'Color', 'g', 'LineStyle', '--', 'linewidth', 1.5);
           set(gca,'YTickLabel',[], 'box', 'off');
        ylabel('');
        xlabel(xlabels(ii));
        xlim(xlimits(ii,:));
        elseif ii<6
            [f,xi] = ksdensity(samples(:,ii));
            plot(xi,f, 'linewidth', 1.5, 'Color', 'black');
            line([mean(samples(:,ii)), mean(samples(:,ii))], [0, max(f)] , 'Color', 'r', 'linewidth', 1.5);
%         line([mean(samples(:,ii)), mean(samples(:,ii))], get(gca, 'ylim') , 'Color', 'r', 'linewidth', 1.5);
        line([true_param(ii), true_param(ii)], [0, max(f)] , 'Color', 'g', 'LineStyle', '--', 'linewidth', 1.5);
        set(gca,'YTickLabel',[], 'box', 'off');
        ylabel('');
        xlabel(xlabels(ii));
        xlim(xlimits(ii,:));
        else
            plot(1,nan, 'black',1,nan, 'r',1,nan, '--g', 'linewidth', 1.5);
            set(hsub, 'Visible', 'off');
            legend("Approx. posterior", "MMSE estimate", "True value");
        end
        
        
       
            
      
        
           
    end
end