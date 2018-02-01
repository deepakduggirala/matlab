function printToFile2(fileId,t,r)
rx = r*1000;
fprintf(fileId, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',t, rx(1,1), rx(2,1), rx(3,1),rx(1,2), rx(2,2), rx(3,2));
end