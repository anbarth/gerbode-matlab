function [ data_matrix] = Apply_collision_wraparound( data_matrix, particle_size, neighbors,right_wall,top_wall)
%takes in a configuration of particles, and if any of them are
%intersecting, resolves the collision by assuming the particles are hard
%disks and places the particles outside of the other's  radius

%data_matrix is the particle positions
%particle_size is the size of the particles in pixels
%neighbors is a list of each particles neighboring particles
%right_wall is the right boundary of the simulation
%top_wall is the top boundary of the simulation




%loop through each aprticle
for i=1:length(neighbors(:,1))
    
    
    
%loops through each neighbor   
    for j=1:length(find(neighbors(i,:)))
        
        
        
        %{
        disp('neighbors is')
        disp(neighbors(i,j))
        disp(length(data_matrix(:,1)))
        %}
        
%get positions of neighboring and current particle        
        r1=data_matrix(neighbors(i,j),[1 2]);
        
        r2=data_matrix(i,[1 2]);
        
%check for wraparound conditions       
        if r1(1) < particle_size && r2(1)>right_wall-particle_size
            r1(1)=r1(1)+right_wall;
        elseif r1(1)>right_wall-particle_size && r2(1) < particle_size
            r1(1)=r1(1)-right_wall;
        end
        
        if r1(2) < particle_size && r2(2)>top_wall-particle_size
            r1(2)=r1(2)+top_wall;
        elseif r1(2) > top_wall - particle_size && r2(2) < particle_size
            r1(2)=r1(2)-top_wall;
        end
        
        
        
        
        %{
        r1=data_matrix(neighbors(i,j),[1 2]);
        r2=data_matrix(i,[1 2]);
        
        %}
 %get distance between particles       
        d=norm(r2-r1,2);
        
        
%if they intersect, perform collison resolution        
        if d <= particle_size
            
            
            
            
            
            
            
            
            
            %This is hard sphere collisions with momentum conversation
            %calculations.  For brownian dynamics, velocity should be 0 and
            %not relevant
            
            
            r=(1/d)*(r2-r1);
            
            v1=data_matrix(neighbors(i,j),[3 4]);
            v2=data_matrix(i,[3 4]);
            
            %switch velocities
            proj1=dot(r, v1)*r;
            pep1=v1-proj1;
            proj2=dot(r, v2)*r;
            pep2=v2-proj2;
            vf1=proj2+pep1;
            vf2=proj1+pep2;
            
            
            data_matrix(neighbors(i,j),[3 4])=vf1;
            data_matrix(i,[3 4])=vf2;
            %repel particles enough so not intersecting
            data_matrix(i,[1 2])=data_matrix(i,[1 2]) + (particle_size-d)/2*r;
            data_matrix(neighbors(i,j),[1 2])=data_matrix(neighbors(i,j),[1 2])-(particle_size-d)/2*r;
            
            
            
            
            
        end
        
    end
end

end

