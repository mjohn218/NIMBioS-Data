AB= dlmread('seed_00001/AB.Cube.dat');
ABC= dlmread('seed_00001/ABC.Cube.dat');
ACB= dlmread('seed_00001/ACB.Cube.dat');
ACCB= dlmread('seed_00001/ACCB.Cube.dat');

complexes = AB+ ABC + ACB + ACCB;

semilogx(complexes(:,1)/4,complexes(:,2));