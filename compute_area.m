function area=compute_area(vertex)
area=0;
    nbe=size(vertex,1)-1;
    if(length(unique(vertex,'rows'))>2 && nbe>2)
        
        % We compute one point inside the cell, to split the cells into triangles
        % to compute the area
        ptK=sum(vertex([1:nbe],:))/nbe;
        
        % Loop over vertices to create the triangles (ptK,vertex1,vertex2)
        for j=1:nbe
            % vertices of the triangle
            v1=vertex(j,:);
            v2=vertex(j+1,:);
            % We check if the cell is indeed star-shaped with respect to this point:
            %	for each edge e with vertices s,s', (ptK-s).n_e should be negative
            %
            %     if ( dot(ptK-v1,(v2-v1)*[0 -1;1 0]) >= 0)
            %         error('Cell %i not star-shaped with respect to point (%f,%f)',i,ptK(1),ptK(2));
            %     end
            %
            % AREA
            % adding area of current triangle
            area = area + 0.5*det([v1-ptK;v2-ptK]);
            
        end
    end
    area=abs(area);


end