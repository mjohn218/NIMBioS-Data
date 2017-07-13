A=0;
R=A;
numfol = 24;

cd /Users/achen70/MCell/AR_Oscillation/react_data

for i = 1:numfol
    
    s = sprintf('%2.2d',i);
    A = A+dlmread(['seed_000' s '/A.Cube.dat']);
    R = R+dlmread(['seed_000' s '/R.Cube.dat']);
    
end

A = A/numfol;
R = R/numfol;

save