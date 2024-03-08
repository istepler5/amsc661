% AMSC 661 Homework 5 1(d) and 1(e)
% Write a code that plots the boundaries of the RAS for the one-parameter
% family of two-set 3rd-order methods for alpha = -1.8,-1.7,...,-1.1 in one
% figure and the boundaries of the RAS for alpha = -1.0,-0.9,...,-0.1 in
% another figure. Include legend in each figure. 
% Will use the boundary locus method

alpha1 = -1.8 : .1 : -1.1;
alpha2 = -1.0 : .1 : -0.1;

RAS_boundary(alpha1,1);
RAS_boundary(alpha2,2);

% part 1(e) - for each value of alpha, check whether the root condition is
% satisfied inside the contours. I chose to evaluate at the single point
% hL = -5+0i for the first set of contours since it is always inside each
% contour, and the point hL=-.25+0i for the second set of contours for the
% same reason. 

where_RAS(alpha1,-5);
where_RAS(alpha2,-.25);


% Function for 1(d)
function [] = RAS_boundary(alpha,fig_num)
    z = linspace(0, 2 * pi, 4000);
    z = exp(z*1i);

    for i = 1:length(alpha)
        rho_coeff = rho(alpha(i));
        sigma_coeff = sigma(alpha(i));

        rho_val = polyval(rho_coeff,z);
        sigma_val = polyval(sigma_coeff,z);

        f = rho_val ./ sigma_val;

        figure(fig_num);
        plot(f,'-','DisplayName',strcat(num2str(alpha(i))));
        hold on;
        xlabel('Re(z)'); ylabel('Im(z)');
        title(['RAS Boundaries for Alpha Set ',num2str(fig_num)]);
    end

    legend('show','Location','northwest');
end

% function for 1(e)
function [] = where_RAS(alpha,hL)
    for i = 1:length(alpha)
        rho_coeff1 = rho(alpha(i));
        sigma_coeff1 = sigma(alpha(i));
    
        % solve for z in rho(z) - hL*sigma(z)=0
        p = rho_coeff1 - hL*sigma_coeff1;
        p_root = roots(p);

        % check the root condition
        abs_p_root = abs(p_root);
        max_root = max(abs_p_root);
        if max_root <= 1
            fprintf('The RAS is inside the contour for alpha0 = %d.\n',...
                alpha(i));
        else
            fprintf('The RAS is outside the contour for alpha0 = %d.\n',...
                alpha(i));
        end
    end
end



function rho_coeff = rho(alpha)
    rho_coeff = [1 alpha -1-alpha];
end

function sigma_coeff = sigma(alpha)
    a = 1/3 - alpha/12;
    b = 4/3 + 2*alpha/3;
    c = 5*alpha/12 + 1/3;
    sigma_coeff = [a b c];
end

