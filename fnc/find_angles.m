function [ AoA AoD ] = find_angles( UE, BS )
% Find AoD and AoA of a user, with respect to a BS
% UE --> User coordinates [x,y]
% BS --> BS coordinates [x,y]

pos_x = -UE(1)+BS(1);
pos_y = -UE(2)+BS(2);

if (pos_x >= 0 && pos_y >=0)
    AoD = atan(pos_y/pos_x);
elseif (pos_x < 0 && pos_y >=0)
    AoD = atan(pos_y/pos_x) + pi;
elseif (pos_x < 0 && pos_y < 0)
    AoD = atan(pos_y/pos_x) + pi;
else
    AoD = atan(pos_y/pos_x) + 2*pi;
end

AoA = -AoD;

end

