function [ image ] = make_dots_wraparound( x0_list,y0_list,r,name,right_wall,top_wall)
%this function draws white particles on a black background
%Inputs:
%x0_list is a list of all the x-coordinates of the particles
%y-_list is a list of all the y-coordiantes of the particles
%r is the radius of the particle in pixels
%name is the name that will be appeneded to any saved images (must contain
%file extension)
%right_wall is where the right most boundary of the frame will be
%top_wall is the top most boundary of where the frame will be
%   Detailed explanation goes here

%initialize matrix of zeros, corresponding to a black background
image=zeros(top_wall,right_wall,3);

cmap = gray(100);
particle_color = cmap(70,:);
edge_color = cmap(40,:);

%add white pixels corresponding to particles
for n=1:length(x0_list)
    %get nth particle position
    x0=x0_list(n);
    y0=y0_list(n);
    
    %define square of which to shade in
    begin_x=ceil(x0-1.1*r);
    end_x=floor(x0+1.1*r);
    begin_y=ceil(y0-1.1*r);
    end_y=floor(y0+1.1*r);
%{   
    if x0<r
        begin_x=1;
    else
        begin_x=ceil(x0-r);
        
    end
    
    if x0>1600 - r
        end_x=1600;
    else
        end_x=floor(x0+r);
        
    end
    
    if y0<r
        begin_y=1;
    else
        begin_y=ceil(y0-r);
    end
    
    if y0>900 - r
        end_y=900;
    else
        end_y=floor(y0+r);
        
    end
    
 %}   
    
    %color in pixels within a distance r of center of particle
    for i= begin_x: end_x
        for j= begin_y : end_y
            if i<1
                x=i+floor(right_wall);
            elseif i>right_wall
                x=i-floor(right_wall);
            else
                x=i;
            end
            
            if j<1
                y=j+floor(top_wall);
            elseif j>top_wall
                y=j-floor(top_wall);
            else
                y=j;
            end
            
          
            
            dist=sqrt((x0 - (i-.5))^2 + (y0 - (j-.5))^2);
           if dist<.9*r
               try
                   image(y,x,:)=particle_color;
               catch
                   disp(x)
                   disp(y)
               end
            elseif dist<r && dist>.9*r
                image(y,x,:)=edge_color;
               
            end
                
            
        end
    end
    
end

%image(435:465,1160:1200)= zeros(31,41);
%ones(31,31,3);

%make an image with matrix
imwrite(image,convertStringsToChars(name));
disp('In make_dots_wraparound-Saved up to:')
disp(name)

end




