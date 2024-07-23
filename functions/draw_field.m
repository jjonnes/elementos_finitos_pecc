function draw_field(nodeCoordinates,elementNodes,...
    elementType,field)
switch elementType
    case 'Q4'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,:),1),nodeCoordinates(elementNodes(i,:),2),field(i,:));
        end
    case 'Q8'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,1:4),1),nodeCoordinates(elementNodes(i,1:4),2),field(i,:));
        end
    case 'Q12'
        for i = 1:size(elementNodes,1)
            patch(nodeCoordinates(elementNodes(i,1:4),1),nodeCoordinates(elementNodes(i,1:4),2),field(i,:));
        end
    otherwise
        disp('Element type not available')
end
end
