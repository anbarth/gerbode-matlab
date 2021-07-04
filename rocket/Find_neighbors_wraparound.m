function [ neighbors_matrix ] = Find_neighbors_wraparound( data_matrix,particle_size, threshold ,right_wall, top_wall)
%takes in data matrix and returns a matrix that is nx20.  Each row
%corresponding to a particle, and the non-zero entries contain the indicies
%of neighbor particles
%Inputs
%data_matrix is the data_matrix with the relevant particle information
%particle_size is the particle+size
%threshold is the multiplier by which we multiply by particle size to see
%if a particle is considered a neighbor or not

 %initialize return variable   
neighbors_matrix=zeros(length(data_matrix(:,1)), 200);

%loop through each particle to find neighbors.  
for i=1:length(data_matrix(:,1))
    
    
    
    x=data_matrix(i,1);
    y=data_matrix(i,2);
    %disp('data_matrix is ')
    %disp(data_matrix)
    if x<threshold*particle_size
        x_neighbors=data_matrix(:,1)<x | data_matrix(:,1)>right_wall - threshold*particle_size;
    elseif x > right_wall-threshold*particle_size
        x_neighbors=data_matrix(:,1)<threshold*particle_size | data_matrix(:,1)>x;
    else 
        x_neighbors=false(length(data_matrix(:,1)),1);
    end
        x_neighbors= x_neighbors | data_matrix(:,1)>x - particle_size*threshold & data_matrix(:,1)< x;
    
    
    if y<threshold*particle_size
        y_neighbors=data_matrix(:,2)<y | data_matrix(:,2)>top_wall - threshold*particle_size;
    elseif y>top_wall-threshold*particle_size
        y_neighbors=data_matrix(:,2)<threshold*particle_size | data_matrix(:,2)>y;
    else 
        y_neighbors=false(length(data_matrix(:,1)),1);
    end
    
        y_neighbors=y_neighbors | data_matrix(:,2)>y - particle_size*threshold & data_matrix(:,2)< y+particle_size*threshold;
    
        
    neighbors=data_matrix(x_neighbors & y_neighbors,7);
    %{
    disp('length is')
    disp(length(neighbors(:,1)))
    %}
    entry= [transpose(neighbors) zeros(1,200-length(neighbors))];
    
    neighbors_matrix(i,:)=entry;

end

