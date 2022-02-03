function domain_plotter(cellpositions,celltype)

[v,c] = voronoin(cellpositions); 

col1 = spring(20);
col2 = parula(18);
col3 = jet(40);
col4 = gray(20);

for ig =1:length(c)
    if celltype(ig)==1;
        if all(c{ig}~=1)% 
        patch(v(c{ig},1),v(c{ig},2),col1(1,:)) % use color i.
            for j = 1:10
            newpoints =v(c{ig},:)+1/(11-j)*(cellpositions(ig,:)-v(c{ig},:));
            patch(newpoints(:,1),newpoints(:,2),col1(j+1+4,:),'EdgeColor','none') % use color i.
            end
        end
  
    elseif celltype(ig)==51;
         if all(c{ig}~=1)% If at least one of the indices is 1, 
          patch(v(c{ig},1),v(c{ig},2),col2(1,:)) % use color i.
            for j = 1:10
            newpoints =v(c{ig},:)+1/(11-j)*(cellpositions(ig,:)-v(c{ig},:));
            patch(newpoints(:,1),newpoints(:,2),col2(j+1,:),'EdgeColor','none') % use color i.
            end
         end
     elseif celltype(ig)==3;
         if all(c{ig}~=1)% If at least one of the indices is 1, 
         patch(v(c{ig},1),v(c{ig},2),col3(1,:)) % use color i.
            for j = 1:10
            newpoints =v(c{ig},:)+1/(11-j+2)*(cellpositions(ig,:)-v(c{ig},:));
            patch(newpoints(:,1),newpoints(:,2),col3(j+29,:),'EdgeColor','none') % use color i.
            end
         end
    elseif celltype(ig)==4;
         if all(c{ig}~=1)% If at least one of the indices is 1, 
         patch(v(c{ig},1),v(c{ig},2),col4(1,:)); % use color i.
            for j = 1:10
            newpoints =v(c{ig},:)+1/(11-j+2)*(cellpositions(ig,:)-v(c{ig},:));
            patch(newpoints(:,1),newpoints(:,2),col4(j+1,:),'EdgeColor','none') % use color i.
            end
         end
         
    end
end
findloc = find(celltype<4);

hold on 
plot(cellpositions(findloc,1),cellpositions(findloc,2),'k.','LineWidth',3)

Healthylocs = find(celltype==4);
MaxlocHealthy_x = max(cellpositions(Healthylocs,1));
MaxlocHealthy_y = max(cellpositions(Healthylocs,2));
xlim([-MaxlocHealthy_x+10,MaxlocHealthy_x-10]);
ylim([-MaxlocHealthy_y+10,MaxlocHealthy_y-10]);

hold off
end