function [ image ] = make_dots_wraparound_color( x0_list,y0_list,r,name,right_wall,top_wall,cmap)
%make_dots_wraparound_color takes an input of x and y positions of
%particles, and draws circular dots at those locations of a specified
%length and color.  it will also use wraparound boundary counditions, so if
%a dot exceeds the screen, the remaining portion will appear on the
%opposite screen.
%Inputs: 
%x0_list is a list of x positions of particles in order
%y0_list is a list of y positions of particles in order
%r is the radius of the particles
%name is the name which the picture will be saved as
%right_wall is the length of the wall to the right where wrapaound boundary
%conditions will be applied
%top_wall is the length of the wall to the top where wraparound boundry
%conditions will be applied
%cmap is a color map, an nx3 matrix with color values to be used to color
%particles


%ADDED BY ELI WEISSLER '19 to fix issue of these images being vertically
%reflected
y0_list = -(y0_list-top_wall/2) + top_wall/2;

image=zeros(top_wall,right_wall,3);

for n=1:length(x0_list)
    x0=x0_list(n);
    y0=y0_list(n);
    
    begin_x=ceil(x0-1.1*r);
    end_x=floor(x0+1.1*r);
    begin_y=ceil(top_wall - (y0+1.1*r));
    end_y=floor(top_wall - (y0-1.1*r));
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
            
          
            
            dist=sqrt((x0 - (i-.5))^2 + (y0 - ((top_wall-j)-.5))^2);
            if dist<r
               try
                   image(y,x,:)=.7*cmap(n,:);
               catch
                   disp(x)
                   disp(y)
               end
            elseif dist>r && dist<1.1*r
                image(y,x,:)=.4*cmap(n,:);
               
            end
                
            
        end
    end
    
end

%image(435:465,460:490,:)=ones(31,31,3);
imwrite(image,convertStringsToChars(name));
disp('In make_dots_wraparound_color-Saved up to')
disp(name)

end




