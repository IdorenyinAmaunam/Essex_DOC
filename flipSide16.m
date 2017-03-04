function flipped = flipSide16(elec16, doflip)

flipped = elec16;
if(doflip)
    flipped = elec16([1 6 5 4 3 2 11 10 9 8 7 16 15 14 13 12]);
end