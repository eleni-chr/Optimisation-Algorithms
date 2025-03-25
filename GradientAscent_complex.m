function GradientAscent_complex
%% NOTE: This function was written using MATLAB R2016a

% NOTE: This function uses a COMPLEX landscape. Another version of this
% function, using a more SIMPLE landscape can be found on GitHub, titled
% "GradientAscent_simple".

% This algorithm contains the function ComplexLandscape which takes an (x,y)
% position as input and returns the height vector of a function of x and y.
% It also contains the function ComplexLandscapeGrad which takes the same
% (x,y) position as input and returns the gradient vector of a function of
% x and y. With respect to the three variables in this algorithm, the
% function will initially:
% 1. use a constant learning rate of 0.1;
% 2. stop after a maximum of 50 iterations of the algorithm;
% 3. as the functions being optimised are known, the exact gradient will be used.
% To perform gradient ascent, the function first calculates the gradient at
% a point and then uses this gradient to update the point it is at. The
% function ComplexLandscapeGrad returns the elements of the gradient of the
% corresponding function as a vector.

%% Plot the landscape
ezmesh(@ComplexLandscape,[-3 7],[-3 7])

% Enter maximum number of iterations of the algorithm, learning rate and mutation range
% Changing these parameters might improve/worsen the performance of the
% algorithm. Try it and see.
NumSteps=50;
LRate=0.1;

% Choose a random starting point with x and y in the range (-3,7) or start from a grid of points (un/comment relevant line)
StartPt=-3+10.*rand(1,2);
%StartPt=point_grid(-3,7,-3,7,1);

% Define how many starting points there are
a=size(StartPt); NumPts=a(1,1); %NumPts is equal to the number of starting points

% Find maximum
GradAscent(StartPt,NumSteps,LRate,NumPts);

end

% Function implementing gradient ascent
function GradAscent(StartPt,NumSteps,LRate,NumPts)
	PauseFlag=1; %PauseFlag is used later to pause to view the output
    heights=zeros(NumSteps,NumPts); %A vector of zeros in which the height of each StartPt will be stored in for later use
	hold on;
	for i=1:NumSteps
        for j=1:NumPts
            % Calculate the 'height' at StartPt
            h=ComplexLandscape(StartPt(j,1),StartPt(j,2));
        
            % Store the height of each StartPt in 'heights' at each iteration for later use
            heights(i,j)=h;
    
            % Plot point(s) on the landscape
            plot3(StartPt(j,1),StartPt(j,2),h,'*','MarkerSize',10)

            % Calculate the gradient at each StartPt
            g=ComplexLandscapeGrad(StartPt(j,1),StartPt(j,2));

            % Calculate the new point and update each StartPt
            StartPt(j,:)=StartPt(j,:) + LRate.*g;

            % Ensure each StartPt is within the specified bounds
            StartPt(j,:)=max([StartPt(j,:);-3 -3]);
            StartPt(j,:)=min([StartPt(j,:);7 7]);
        end

        % Pause to view output
	    if(PauseFlag)
	        x=input('Press return to continue\nor 0 and return to stop pausing\n');
	        if(x==0)
                PauseFlag=0;
            end
        end
    end
	hold off

    %Plot a pseudocolour plot
    if NumPts>1
        m=max(heights); %finds the maximum value of h reached by each of the starting points
        Colours=vec2mat(m,sqrt(NumPts)); %The values of the elements of Colours specify the color in each cell of figure 2
        figure(2) %a pseudocolour plot with the darkest cells indicating points that reached a larger height
        pcolor(Colours); axis off; colormap(flipud(bone)); shading faceted; h=colorbar; set(h,'FontSize',12); title('Gradient Ascent','FontSize',12);
    end
    max_height_reached=max(max(heights)) %returns the maximum height that was reached by the algorithm
end

% Definition of Complex landscape
function [f]=ComplexLandscape(x,y)
	f=4*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*(x/5 - x^3 - y^5)*exp(-x^2-y^2) -(1/3)*exp(-(x+1)^2 - y^2)-1*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9)*exp(-(x-3)^2-(y-3)^2);
end

% Definition of gradient of Complex landscape
function [g]=ComplexLandscapeGrad(x,y)
	g(1)=-8*exp(-(x^2)-(y+1)^2)*((1-x)+x*(1-x)^2)-15*exp(-x^2-y^2)*((0.2-3*x^2) -2*x*(x/5 - x^3 - y^5)) +(2/3)*(x+1)*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*(14*(x-3)^6-2*(x-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
	g(2)=-8*(y+1)*(1-x)^2*exp(-(x^2)-(y+1)^2) -15*exp(-x^2-y^2)*(-5*y^4 -2*y*(x/5 - x^3 - y^5)) +(2/3)*y*exp(-(x+1)^2 - y^2)-1*exp(-(x-3)^2-(y-3)^2)*((-1.5*(y-4)^4+9*(y-3)^8)-2*(y-3)*(2*(x-3)^7-0.3*(y-4)^5+(y-3)^9));
end

%Function creating a grid of points on the landscape
function M = point_grid(minx,maxx,miny,maxy,step)
x=minx;
y=miny;
M=zeros(((maxx-miny)./step).^2,2);
i=1;
while x<=maxx
    y=miny;
    while y<=maxy
        M(i,:)=[x y];
        i=i+1;
        y=y+step;
    end
    x=x+step;
end
end