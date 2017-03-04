function flipped = flipSide64(elec64, doflip)

flipped = elec64;
if(doflip)
    flipped = elec64([3:-1:1 7 6 5 4 62 16:-1:8 61 25:-1:17 34:-1:26 43:-1:35 52:-1:44 57:-1:53 60:-1:58]);
end